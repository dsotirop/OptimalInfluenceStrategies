% This script file performs all necessary symbolic computations involved 
% in order to address the problem of optimal influence for network with one
% company and three consumers. The social interaction matrix is given by T.


% Clear workspace and command window.
clc
clear all

% Setup symbolic variables. 
syms T1 T2 T3 P1 P2 P3 
syms theta1 theta2 theta3
syms VS0 VS1 VS2 VS3

% Symbolic variables interpretation.
% T1: T_10
% T2: T_20
% T3: T_30
% P1: P1(0)
% P2: P2(0)
% P3: P3(0)
% theta1: T_21 / T_12
% theta2: T_31 / T_13
% theta3: T_32 / T_23


% Setup constant values for corresponding symbolic variables.
P1_const = 0.1;
P2_const = 0.2;
P3_const = 0.3;
theta1_const = 0.1;
theta2_const = 0.1;
theta3_const = 0.1;
Co = 0.2;
% Set minimum and maximim values for the variables of the optimization
% problem.
T1_min = 0;
T2_min = 0;
T3_min = 0;
T1_max = 1/3;
T2_max = (2-theta1_const) / 3;
T3_max = (3-theta2_const-theta3_const) / 3;

% Report Experimental Setup.
fprintf('-------------------------------------------------------------\n');
fprintf('Experimental Setup:\n');
fprintf('P1(0) = %d\n',P1_const);
fprintf('P2(0) = %d\n',P2_const);
fprintf('P3(0) = %d\n',P3_const);
fprintf('Theta1 = %d\n',theta1_const);
fprintf('Theta2 = %d\n',theta2_const);
fprintf('Theta3 = %d\n',theta3_const);
fprintf('Co = %d\n',Co);
fprintf('-------------------------------------------------------------\n');
fprintf('Constraints on T1, T2 and T3:\n');
fprintf('0 < T1 < %f\n',T1_max);
fprintf('0 < Ô2 < %f\n',T2_max);
fprintf('0 < Ô3 < %f\n',T3_max);

% Assign constraints on symbolic variables T1, T2 and T3.
assume(and((T1>0),(T1<T1_max)));
assume(and((T2>0),(T2<T2_max)));
assume(and((T3>0),(T3<T3_max)));

% Setup social interaction matrix components.
T = [1/2 1/6 1/6 1/6;...
     T1 (1/3)-T1 1/3 1/3;...
     T2 theta1/3 ((2-theta1)/3)-T2 1/3;...
     T3 theta2/3 theta3/3 ((3-theta2-theta3)/3)-T3;];
 
% Setup the components of the final influence vector.
VS = [VS0 VS1 VS2 VS3];
% Setup the corresponding system initiating from the
% left eigenvalues computation problem.
Y = VS * (T - eye(4));
System1 = [Y(1)==0;Y(2)==0;Y(3)==0;VS0+VS1+VS2+VS3==1];
% Solve the corresponding linear system.
Solution1 = solve(System1,VS0,VS1,VS2,VS3);
S0 = Solution1.VS0;
S1 = Solution1.VS1;
S2 = Solution1.VS2;
S3 = Solution1.VS3;
S = [S0;S1;S2;S3];
% Compute the limiting influence vector components.
B = [1 P1 P2 P3];
X = B * S;
X = simplify(collect(expand(X)));
% Compute corresponding partial derivatives with respect to T1, T2 and T3.
F1 = diff(X,T1);
F2 = diff(X,T2);
F3 = diff(X,T3);
F1 = simplify(collect(expand(F1)));
F2 = simplify(collect(expand(F2)));
F3 = simplify(collect(expand(F3)));

% Substitute some of the symbolic variables with corresponding constant
% values.
G1 = subs(F1,[P1 P2 P3 theta1 theta2 theta3],[P1_const,P2_const,P3_const,theta1_const,theta2_const,theta3_const]);
G2 = subs(F2,[P1 P2 P3 theta1 theta2 theta3],[P1_const,P2_const,P3_const,theta1_const,theta2_const,theta3_const]);
G3 = subs(F3,[P1 P2 P3 theta1 theta2 theta3],[P1_const,P2_const,P3_const,theta1_const,theta2_const,theta3_const]);
G1 = simplify(collect(expand(G1)));
G2 = simplify(collect(expand(G2)));
G3 = simplify(collect(expand(G3)));

% Solve the corresponding nonlinear system utilizing a symbolic solver.
System2 = [G1==Co;G2==Co;G3==Co];
Solution2 = solve(System2,T1,T2,T3);

% Generate the function m file that will be subsequently utilized in the
% numeric optimization problem.
GenerateNonLinearSystemFunction('NonLinearSystem',System2,'T');

% Solve the corresponding nonlinear system utilizing a numeric solver.
EquationsNum = length(System2);
N = 20;
lb = [T1_min,T2_min,T3_min];
ub = [T1_max,T2_max,T3_max];
Solutions = SolveNonLinearSystemFmincon(EquationsNum,N,lb,ub);