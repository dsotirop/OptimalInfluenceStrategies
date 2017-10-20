function [Solutions,ExitFlags] = SolveNonLinearSystem(EquationsNum,N)

% This function performs numerical optimization on the nonlinear
% optimization problem defined in the function file "NonLinearSystem"
% by utilizing the fmincon routine.

rng default 
pts = 100*rand(N,EquationsNum);
soln = zeros(N,EquationsNum);
exflags = zeros(N,EquationsNum);
opts = optimoptions(@fsolve,'Display','off','Algorithm','trust-region-reflective');
for k = 1:N
    [soln(k,:),~,exflags(k,:)] = fsolve(@NonLinearSystem,pts(k,:),opts);
end;
Solutions = soln;
ExitFlags = exflags;
end