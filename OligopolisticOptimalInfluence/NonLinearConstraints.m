function [c,ceq] = NonLinearConstraints(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% This function defines the sets of nonlinear equality and inequality
% constraints that are to be associated with the problem of determining the 
% optimal influence strategies within the context of oligopolistic markets.

% For this particular problem there will be no nonlinear inequality
% constraints.
c = [];

% The set of nonlinear equality constraints, however, will be determined by
% the code function "NonLinearEquationsSystem.m" which defines the system
% of nonlinear equations that have to be satisfied.
ceq = NonLinearEquationsSystem(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);


end

