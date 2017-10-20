% This script file performs fundamental symbolic operations in order to
% assist the process of obtaining an analytical solution for the problem of
% optimal influences estimation on a network comprising of one company and
% two consumers. In particular, this piece of code focuses on the problem
% arising when T1+T2 = 1.

clc
clear all

% Define the underlying symbolic variables.
syms T1 T2 Lambda1 Lambda2 Theta1 Theta2 A C Gamma P1 P2

% Define the fundamental problem quantities.
Q = T1*T2 + (Lambda2+Theta1)*T1 + (Lambda1+Theta2)*T2 + (Lambda1+Lambda2)*(Theta1+Theta2);
S0 = (T1*T2 + T1*Theta1 + T2*Theta2) / Q;
S1 = (Lambda1*T2 + (Lambda1+Lambda2)*Theta1) / Q;
S2 = (Lambda2*T1 + (Lambda1+Lambda2)*Theta2)/ Q;
S = [S0 S1 S2];
P = [1 P1 P2];
X = S * P.';
Fo = Gamma * (T1+T2) - (1/2)*(A-C)^2*X;

% In order to perform symbolic computations in a more readable form the
% following quantities may be defined.
F0 = T1*T2 + T1*Theta1 + T2*Theta2;
F1 = P1 * (Lambda1*T2 + (Lambda1+Lambda2)*Theta1);
F2 = P2 * (Lambda2*T1 + (Lambda1+Lambda2)*Theta2);
F = F0 + F1 + F2;
% Given that T2 = 1 - T1 we get:
F = subs(F,T2,1-T1);
Q = subs(Q,T2,1-T1);
% Get the coefficients of the polynomials.
CF = coeffs(F,T1);
CQ = coeffs(Q,T1);

A2 = CF(3);
A1 = CF(2);
A0 = CF(1);

B2 = CQ(3);
B1 = CQ(2);
B0 = CQ(1);

% Compute the G(T) polynomial.
syms a2 a1 a0
syms b2 b1 b0
syms T

L = (2*a2*T+a1)*(b2*T^2+b1*T+b0);
R = (2*b2*T+b1)*(a2*T^2+a1*T+a0);
L = collect(L,T);
R = collect(R,T);
G = L - R;
G = collect(G,T);
CG = coeffs(G,T);

Ga = CG(3);
Gb = CG(2);
Gc = CG(1);

Ga = subs(Ga,[a2 b2],[A2 B2]);
Gb = subs(Gb,[a2 b2],[A2 B2]);
Gc = subs(Gc,[a2 b2],[A2 B2]);

Ga = subs(Ga,[a1 a0 b1 b0],[A1 A0 B1 B0]);
Gb = subs(Gb,[a1 a0 b1 b0],[A1 A0 B1 B0]);
Gc = subs(Gc,[a1 a0 b1 b0],[A1 A0 B1 B0]);

DG = Gb^2 - 4*Ga*Gc;
T_a =(-Gb+sqrt(DG))/(2*Ga);
T_b =(-Gb-sqrt(DG))/(2*Ga);

% When A1=B1 (that is Lambda1*(P1-1) = Lambda2*(P2-1)
% then G(T) will not be a second degree polynomial and therefore
% T1_opt will be given by:
%            A0*B1 - A1*B0 
% T1_opt =  -------------- 
%             2*(A0-B0)

syms Lambda P Theta K

T1_optimal = (A0*B1 - A1*B0)/(2*(A0-B0));
T1_optimal = subs(T1_optimal,[Lambda1 Lambda2 P1 P2],[Lambda Lambda P P]);
T1_optimal = collect(T1_optimal,Theta1);

T1_optimal = (Theta1 - Theta2 + K)/2;

S0_optimal = subs(S0,[T2 Lambda1 Lambda2 P1 P2],[K-T1,Lambda,Lambda,P,P]);
S0_optimal = subs(S0_optimal,T1,T1_optimal);

S1_optimal = subs(S1,[T2 Lambda1 Lambda2 P1 P2],[K-T1,Lambda,Lambda,P,P]);
S1_optimal = subs(S1_optimal,T1,T1_optimal);

S2_optimal = subs(S2,[T2 Lambda1 Lambda2 P1 P2],[K-T1,Lambda,Lambda,P,P]);
S2_optimal = subs(S2_optimal,T1,T1_optimal);

X_optimal = subs(X,[T2 Lambda1 Lambda2 P1 P2],[K-T1,Lambda,Lambda,P,P]);
X_optimal = subs(X_optimal,T1,T1_optimal);

Fo_optimal = subs(Fo,[T2 Lambda1 Lambda2 P1 P2],[K-T1,Lambda,Lambda,P,P]);
Fo_optimal = subs(Fo_optimal,T1,T1_optimal);

DS0_optimal1 = diff(S0_optimal,Theta1);
DS0_optimal2 = diff(S0_optimal,Theta2);

DS1_optimal_theta1 = diff(S1_optimal,Theta1);
DS1_optimal_theta2 = diff(S1_optimal,Theta2);
DS1_optimal_lambda = diff(S1_optimal,Lambda);

DX_optimal_theta1 = diff(X_optimal,Theta1);
DX_optimal_theta2 = diff(X_optimal,Theta2);
DX_optimal_lambda = diff(X_optimal,Lambda);

DS0_optimal_K = diff(S0_optimal,K);
DS1_optimal_K = diff(S1_optimal,K);
DS2_optimal_K = diff(S2_optimal,K);
DX_optimal_K = diff(X_optimal,K);

DS0_optimal_K = factor(DS0_optimal_K);
DS1_optimal_K = factor(DS1_optimal_K);
DS2_optimal_K = factor(DS2_optimal_K);
DX_optimal_K = factor(DX_optimal_K);

LS0_optimal_K = latex(DS0_optimal_K);
LS1_optimal_K = latex(DS1_optimal_K);
LS2_optimal_K = latex(DS2_optimal_K);
LX_optimal_K = latex(DX_optimal_K);

DFo_optimal_K = diff(Fo_optimal,K);

% When A1-B1<>0 and A0-B0<>0 (but Lambda1=Lambda2=Lambda)
GA = subs(Ga,[Lambda1 Lambda2],[Lambda Lambda]);
GB = subs(Gb,[Lambda1 Lambda2],[Lambda Lambda]);
GC = subs(Gc,[Lambda1 Lambda2],[Lambda Lambda]);
GD = GB^2 - 4*GA*GC;
Ta = (-GB + sqrt(GD))/(2*GA);
Tb = (-GB - sqrt(GD))/(2*GA);