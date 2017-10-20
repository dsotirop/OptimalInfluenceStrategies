function [F_A,F_B,F_A_Rev,F_A_Cost,F_B_Rev,F_B_Cost] = OligopolisticFOptimal(T1_A,T2_A,T1_B,T2_B,pA,pB,Q_A,Q_B,C,Gamma)

% This function computes the optimal profits (F_A and F_B) for firms A and
% B with respect to the optimal:
% i)   connection strengths T1_A, T2_A, T1_B and T2_B 
% ii)  prices pA and pB 
% iii) the optimal quantities Q_A and Q_B 
% iv)  the paremeters C and Gamma.

% Get the dimensionality of the input vectors.
Lo = length(T1_A);

if (Lo==1)
    % Compute the optimal value for the profit of firm A.
    F_A =  Q_A * (pA - C) - Gamma*(T1_A^2 + T2_A^2);
    F_A_Rev = Q_A * (pA - C);
    F_A_Cost = Gamma*(T1_A^2 + T2_A^2);
    % Compute the optimal value for the profit of firm B.
    F_B = Q_B * (pB - C) - Gamma*(T1_B^2 + T2_B^2);
    F_B_Rev = Q_B * (pB - C);
    F_B_Cost = Gamma*(T1_B^2 + T2_B^2);
else
    % Compute the optimal value for the profit of firm A.
    F_A =  Q_A .* (pA - C) - Gamma.*(T1_A.^2 + T2_A.^2);
    F_A_Rev = Q_A .* (pA - C);
    F_A_Cost = Gamma.*(T1_A.^2 + T2_A.^2);
    % Compute the optimal value for the profit of firm B.
    F_B = Q_B .* (pB - C) - Gamma.*(T1_B.^2 + T2_B.^2);
    F_B_Rev = Q_B .* (pB - C);
    F_B_Cost = Gamma.*(T1_B.^2 + T2_B.^2);
end;

end

