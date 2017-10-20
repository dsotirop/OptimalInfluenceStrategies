% This script file performs fundamental symbolic operations in order to
% assist the process of obtaining an analytical solution for the problem of
% optimal influences estimation on a network comprising of one company and
% two consumers. In particular, this piece of code focuses on the problem
% arising when T1+T2 <= c.

clc
clear all

% Define fundamental symbolic variables.
syms T1 T2 A C Gamma P1 P2 Lambda1 Lambda2 Theta1 Theta2

% Define fundamental problem quantities.
F = T1*T2 + (Lambda2*P2+Theta1)*T1 + (Lambda1*P1+Theta2)*T2 + (Lambda1+Lambda2)*(Theta1*P1+Theta2*P2);
Q = T1*T2 + (Lambda2+Theta1)*T1 +(Lambda1+Theta2)*T2 + (Lambda1+Lambda2)*(Theta1+Theta2);
Fo = Gamma*(T1+T2) - ((A-C)^2*(F/Q))/2;

DF_T1 = diff(F,T1);
DF_T2 = diff(F,T2);
DQ_T1 = diff(Q,T1);
DQ_T2 = diff(Q,T2);

G1 = DF_T1*Q - F*DQ_T1; % Proves to be a second degree polynomial of T2.
G2 = DF_T2*Q - F*DQ_T2; % Proves to be a second degree polynomial of T1.

G1 = simplify(collect(expand((G1)))); % G1 = A2 * T2^2 + A1*T2 + A0.
G2 = simplify(collect(expand((G2)))); % G2 = B2 * T1^2 + B1*T1 + B0.

CG1 = coeffs(G1,T2);
CG2 = coeffs(G2,T1);

A2 = CG1(3);
A1 = CG1(2);
A0 = CG1(1);

B2 = CG2(3);
B1 = CG2(2);
B0 = CG2(1);

% Consider the case where Lambda1 = Lambda2 = Lambda.
syms Lambda
A2 = subs(A2,[Lambda1 Lambda2],[Lambda Lambda]);
A1 = subs(A1,[Lambda1 Lambda2],[Lambda Lambda]);
A0 = subs(A0,[Lambda1 Lambda2],[Lambda Lambda]);

B2 = subs(B2,[Lambda1 Lambda2],[Lambda Lambda]);
B1 = subs(B1,[Lambda1 Lambda2],[Lambda Lambda]);
B0 = subs(B0,[Lambda1 Lambda2],[Lambda Lambda]);

syms x y
g1 = A2*y^2 + A1*y + A0;
g2 = B2*x^2 + B1*x + B0;
G = g1 - g2;
Sol1 = solve(G==0,x,y);
Sol1.x
Sol1.y

Qo = x*y + (Lambda2+Theta1)*y +(Lambda1+Theta2)*x + (Lambda1+Lambda2)*(Theta1+Theta2);
Qo = subs(Qo,[Lambda1 Lambda2],[Lambda Lambda]);
syms z W
g1 = subs(g1,y,Sol1.y(1));
g2 = subs(g2,x,Sol1.x(1));
Qo = subs(Qo,[x y],[Sol1.x(1),Sol1.y(1)]);

Sol2 = solve(g1==W*Qo^2,z);