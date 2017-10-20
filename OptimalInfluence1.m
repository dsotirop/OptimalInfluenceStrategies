% This script file performs all necessary symbolic computations involved 
% in order to address the problem of optimal influence for network with one
% company and two consumers. The social interaction matrix is given by T.

% Clear workspace and command window.
clc
clear all

% Setup symbolic variables. 
syms T1 T2 P1 P2 theta
syms VS0 VS1 VS2
syms C A P G

% Symbolic variables interpretation.
% T1: T_10
% T2: T_20
% P1: P1(0)
% P2: P2(0)
% theta1: T_21 / T_12

% Setup constant values for corresponding symbolic variables.
P1_const = 0.2;
P2_const = 0.2; 
theta_const = 0.0;
%Co = 0.2;
A_const = 3;
C_const = 2;
G_const = 0.1;
P_const = (A_const + C_const)/2;
Co = (2*G_const)/((A_const-C_const)^2);

% Set minimum and maximim values for the variables of the optimization
% problem.
T1_min = 0;
T2_min = 0;
T1_max = 1/2;
T2_max = 1-(theta_const/2);

% Report Experimental Setup.
fprintf('-------------------------------------------------------------\n');
fprintf('Experimental Setup:\n');
fprintf('P1(0) = %d\n',P1_const);
fprintf('P2(0) = %d\n',P2_const);
fprintf('Theta = %d\n',theta_const);
fprintf('Co = %d\n',Co);
fprintf('-------------------------------------------------------------\n');
fprintf('Constraints on T1 and T2:\n');
fprintf('0 < T1 < %f\n',T1_max);
fprintf('0 < Ô2 < %f\n',T2_max);

% Assign constraints on symbolic variables T1, T2 and T3.
assume(and((T1>=0),(T1<=T1_max)));
assume(and((T2>=0),(T2<=T2_max)));

% Setup social interaction matrix components.
T = [1/2 1/4 1/4;...
     T1 ((1/2) - T1) 1/2;...
     T2 theta/2 1-(theta/2)-T2];

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
B = [1 P1 P2];
X = B * S;
X = simplify(collect(expand(X)));

% Setup the objective function to be miminimized.
Fo = -2*(P-C)*X*(A-P) + G*(T1+T2);
HFo = hessian(Fo,[P,T1,T2]);
Go = subs(Fo,[P1 P2 theta A P C G],[P1_const,P2_const,theta_const,A_const,P_const,C_const,G_const]);
C_0 = 2*G_const / ((A_const-C_const)^2);
C_0_string = num2str(C_0);
P1_string = num2str(P1_const);
P2_string = num2str(P2_const);
theta_string = num2str(theta_const);
name_string = strcat(['P1(0) = ' P1_string ' P2(0) = ' P2_string ' theta = ' theta_string ' Co = ' C_0_string]);
figure('Name',name_string)
ezsurf(Go,[0,1,0,1])
colormap(jet)
legend('Go');

% Setup the first component of the objective function to be minimized.
F = (((A+C)/2)-C) * X * (A-((A+C)/2));
F = simplify(collect(F));
% Compute the corresponding Hessian matrix.
HF = hessian(F,[T1 T2]);

% Compute corresponding partial derivatives with respect to T1 and T2.
F1 = diff(X,T1);
F2 = diff(X,T2);
F1 = simplify(collect(expand(F1)));
F2 = simplify(collect(expand(F2)));

% Substitute some of the symbolic variables with corresponding constant
% values.
G1 = subs(F1,[P1 P2 theta],[P1_const,P2_const,theta_const]);
G2 = subs(F2,[P1 P2 theta],[P1_const,P2_const,theta_const]);
G1 = simplify(collect(expand(G1)));
G2 = simplify(collect(expand(G2)));

% Solve the corresponding nonlinear system utilizing a symbolic solver.
System2 = [G1==Co;G2==Co];
Solution2 = solve(System2,T1,T2);

% Generate the function m file that will be subsequently utilized in the
% numeric optimization problem.
GenerateNonLinearSystemFunction('NonLinearSystem',System2,'T');

% Solve the corresponding nonlinear system utilizing a numeric solver.
EquationsNum = length(System2);
N = 10;
lb = [T1_min,T2_min];
ub = [T1_max,T2_max];
[Solutions,Fvals,ExitFlags] = SolveNonLinearSystemFmincon(EquationsNum,5*N,lb,ub);

% Solve the corresponding unconstrained optimization problem.
[UnconsSolutions,UnconsExitFlags] = SolveNonLinearSystem(EquationsNum,5*N);

% Generate a graphical representation of the unconstrained optimization
% problem to be solved.
H1 = G1 - Co;
H2 = G2 - Co;
figure()
ezsurf(H1,[0,1,0,1])
colormap(jet)
freezeColors
hold on
ezsurf(H2,[0,1,0,1])
colormap(copper)
legend('H1','H2');
hold off

NumericalSolution = Solutions(1,:);
T1_opt = NumericalSolution(1);
T2_opt = NumericalSolution(2);

S0_opt = subs(S0,[T1 T2 theta],[T1_opt T2_opt theta_const]);
S1_opt = subs(S1,[T1 T2 theta],[T1_opt T2_opt theta_const]);
S2_opt = subs(S2,[T1 T2 theta],[T1_opt T2_opt theta_const]);

S_opt = eval([S0_opt S1_opt S2_opt])

% Checj results with newer code version.
Delta_const = A_const - C_const;
[T1_opt,T2_opt,Fval] = OptimalInfluences(P1_const,P2_const,theta_const,Delta_const,G_const)