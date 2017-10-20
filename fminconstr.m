function [c,ceq] = fminconstr(x)

c = []; % no nonlinear inequality
ceq = NonLinearSystem(x); % the fsolve objective is fmincon constraints