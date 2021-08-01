% This script file attempts to investigate the conditions are which the
% first stage optimal quantity is positive with respect to the underlying
% parameters M, K, C, XA and XB.

% We have shown that the first stage optimal quantities are given by the
% following equations:
% QA_opt = alpha * XA + beta * XB - gamma' [1]
% QB_opt = beta * XA + alpha * XB - gamma' [2]
% where the auxiliary parameters alpha, beta, gamma and gamma' are given by
%
%           K*M - 2                2*K - M                    C  
% alpha = ----------- [3] beta = ----------- [4] gamma = ----------- [5]
%           M^2 - 4                M^2 - 4                  M - 2
%
% gamma' = gamma * (M - 1)[6]
%

clc
clear all

% Define the necessary symbolic variables.
syms XA XB alpha beta gamma_prime K M C

% Add fundamental range constraints on the symbolic variables.
assume(C>0);
assume(and((0<K),(K<1)));
assume(and((0<M),(M<1)));
assume(and((0<XA),(XA<1)));
assume(and((0<XB),(XB<1)));

alpha = (K*M-2)/(M^2-4);
beta = (2*K-M)/(M^2-4);
gamma_prime = (C*(M-1))/(M-2);

% Define the expressions for the first stage optimal quantiies.
QA_opt = alpha * XA + beta * XB - gamma_prime;
QB_opt = beta * XA + alpha * XB - gamma_prime;

% Obtain the fractional components NQA, DQA and NQB, DQB for the two
% expressions by computing the corresponding numerators and denumerators.
[NQA,DQA] = numden(QA_opt);
[NQB,DQB] = numden(QB_opt);

% Collect the terms and coefficients for the numerators of the expressions
% QA_opt and QB_opt with respect to M.
[TNQA,CNQA] = coeffs(NQA,M);
[TNQB,CNQB] = coeffs(NQB,M);

% Define the disriminants of the expressions NQA and NQB.
DNQA = TNQA(2)^2 - 4*TNQA(1)*TNQA(3);
DNQB = TNQB(2)^2 - 4*TNQB(1)*TNQB(3);

% Rexpress the discriminant-related quantities as polynomials of C.
[TDNQA,CDNQA] = coeffs(DNQA,C);
[TDNQB,CDNQB] = coeffs(DNQB,C);

% According to the previous derivations the numerators of the quantities QA
% and Qb (first-stage optimality) may be written by considering the following
% symbols W = KXA - XB and Z = KXB - XA
syms W Z
% as:
GA = -C*M^2 + (W-C)*M + 2*(Z+C);
GB = -C*M^2 + (Z-C)*M + 2*(W+C);

% Therefore, the determinants of GA and GB with respect to M will be given
% with the following code.
[TGA,CGA] = coeffs(GA,M);
[TGB,CGB] = coeffs(GB,M);
DGA = TGA(2)^2 - 4*TGA(1)*TGA(3);
DGB = TGB(2)^2 - 4*TGB(1)*TGB(3);
% Express discriminants DGA and DGB as polynomials of C.
[TDGA,CDGA] = coeffs(DGA,C);
[TDGB,CDGB] = coeffs(DGB,C);
% Quantities DGA and DGB as functions of the parameter C will be denoted as
% FA and FB such that:
% FA(C) = 9C^2 + (8Z-2W)C + W^2
% FB(C) = 9C^2 + (8W-2Z)C + W^2
% Compute the discriminatns of polynomials FA and FB.
DFA = TDGA(2)^2 - 4*TDGA(1)*TDGA(3);
DFB = TDGB(2)^2 - 4*TDGB(1)*TDGB(3);

% Express discriminants DFA and DFB as functions of the parameter W and Z 
% respectively.
[TDFA,CDFA] = coeffs(DFA,W);
[TDFB,CDFB] = coeffs(DFB,Z);

% Quantities DFA and DFB as functions of the parameters W and Z
% respectively will be denoted as:
% HA(W) = -32W^2 - 32WZ + 64Z^2
% HB(Z) = -32Z^2 - 32ZW + 64W^2
% Compute the discriminants of the polynomials HA and HB.
DHA = TDFA(2)^2 - 4*TDFA(1)*TDFA(3);
DHB = TDFB(2)^2 - 4*TDFB(1)*TDFB(3);