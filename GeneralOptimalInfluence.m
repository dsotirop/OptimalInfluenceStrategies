% This script file performs all necessary symbolic computations involved 
% in order to address the problem of optimal influence for network with one
% company and two consumers. The social interaction matrix is given by T.
% This time, however, the social interaction matrix T is in its most
% general form for a network with 3 agents (1 enterprise and 2 consumers).

% Clear workspace and command window.
clc
clear all

% Setup symbolic variables.
syms Lambda1 Lambda2
syms T1 T2
syms Theta1 Theta2
syms P0 P1 P2
syms VS0 VS1 VS2
syms Delta Gamma

% Symbolic variables interpretation.
% Notice that parameters Theta1 and Theta2 do not correspond to relative
% influences according to the current formuluation.
% T1: T_10
% T2: T_20
% P0: P0(0)
% P1: P1(0)
% P2: P2(0)
% Lambda1: T_01
% Lambda2: T_02
% Theta1: T_21
% Theta2: T_12
% Delta: Corresponds to the difference (A - C)
% Gamma: Corresponds to the Gamma parameter of the model.

% Setup social interaction matrix components.
T = [1-Lambda1-Lambda2 Lambda1 Lambda2;...
     T1 1-T1-Theta2 Theta2;...
     T2 Theta1 1-T2-Theta1];
% Setup the components of the final influence vector. 
VS = [VS0 VS1 VS2];
% Setup the corresponding system initiating from the
% left eigenvalues computation problem.
Y = VS * (T - eye(3));
System1 = [Y(1)==0;Y(2)==0;VS0+VS1+VS2==1];
% Solve the corresponding linear system.
Solution1 = solve(System1,VS0,VS1,VS2);
S0 = Solution1.VS0;
S1 = Solution1.VS1;
S2 = Solution1.VS2;
S = [S0;S1;S2];
% Compute the limiting influence vector components.
B = [P0 P1 P2];
X = B * S;
X = simplify(collect(expand(X)));
% Check that the sum of S components equals 1.
Ssum = simplify(collect(expand(sum(S))));
if(Ssum==1)
    fprintf('Elements of S sum up to 1\n');
end;
% Setup the objective function to be miminimized.
Fo = -(1/2)*Delta^2*X + Gamma*(T1+T2);