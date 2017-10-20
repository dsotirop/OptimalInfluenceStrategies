% This script file sets up an optimal ifluence problem which is solved as
% a non-linear minimization problem.

% Set constant values for the optimization problem.
P1_const = 0.4;
P2_const = 0.4; 
theta_const = 0.8;
A_const = 3;
C_const = 1;
G_const = 1;
P_const = (A_const + C_const)/2;
Co = (2*G_const)/((A_const-C_const)^2);

% Set boundary values for the free parameters of the optimization
% problem.
T1_min = 0;
T1_max = 1/2;
T2_min = 0;
T2_max = 1 - (theta_const/2);

% Report Experimental Setup.
fprintf('-------------------------------------------------------------\n');
fprintf('Experimental Setup:\n');
fprintf('P1(0) = %d\n',P1_const);
fprintf('P2(0) = %d\n',P2_const);
fprintf('Theta = %d\n',theta_const);
fprintf('A = %d\n',A_const);
fprintf('C = %d\n',C_const);
fprintf('G = %d\n',G_const);
fprintf('Co = %d\n',Co);
fprintf('-------------------------------------------------------------\n');
fprintf('Constraints on T1 and T2:\n');
fprintf('0 < T1 < %f\n',T1_max);
fprintf('0 < Ô2 < %f\n',T2_max);

% Set lower and upper bounds for the free parameters of the optimization.
lb = [T1_min,T2_min];
ub = [T1_max,T2_max];

% Find minimum of the corresponding nonlinear function for a set of 
% N different initial points.
N = 20;
Dimensionality = 2;
[Solutions,Fvals,ExitFlags] = MinimizeNonLinearFunctionFmincon(Dimensionality,N,lb,ub,P1_const,P2_const,theta_const,A_const,C_const,G_const);
