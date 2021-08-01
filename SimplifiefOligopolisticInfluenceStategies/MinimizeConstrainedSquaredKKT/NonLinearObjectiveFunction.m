function [F] = NonLinearObjectiveFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function defines the most general form of the system of nonlinear
% equations that has to be solved in order to obtain the optimal influence
% strategies within the simplified oligopolostic environment of the two
% firms {Firm A, Firm B} and one consumer {C}. However, the system of
% nonlinear equations will be reformulated as a nonlinear function
% optimization problem.

% The input vector T is assumed to be of the form:
% T = [TA TB]

% Get the influence related variables (the actual optimization varibles).
TA = T(1);
TB = T(2);

% Mind that FIRST or SECOND ALTERNATIVE code blocks should be followed.
% -------------------------------------------------------------------------
%                          FIRST ALTERNATIVE
% -------------------------------------------------------------------------
% The system of nonlinear equations corresponds to the first order
% conditions on the profit functions fa and fb for products A and B. In 
% fact, we have a nonlinear system of the following form:
%          [F1 F2] = [F(1) F(2)] =  [0 0] [I]
% where
%               d fa                  d fb
%       F(1) = --------  and F(2) = -------- [II] 
%               d TA                  d Tb
%F(1) = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
%F(2) = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
% -------------------------------------------------------------------------
%                          SECOND ALTERNATIVE
% -------------------------------------------------------------------------
% Since the first derivatives of fa(TA,TB) and fb(TA,TB) may be expressed
% as fractions of the form:
%          d fa      Ua(TA,TB)              d fb       Ub(TA,TB)  
% F(1) = -------- = ----------  and F(2) = --------  = --------- [III] 
%          d TA      Va(TA,TB)              d Tb       Vb(TA,TB)   
%
% where Va(TA,TB) = Vb(TA,TB) <> 0, the system of nonlinear equations may
% be reduced to a system of (4th degree) polynomial equations with respect
% to TA and TB. Actually, Ua is a fourth degree polynomial with respect to
% TA where the corresponding coefficients are parameterized by TB, Caa(TB).
% This is also true for Ub whicg is a fourth degree polynomial with respect
% to TB where the corresponding coefficients are parameterized by TA,
% Cbb(TA). Thus, the system of polynomial equaitions will be given by:
%                       (i):  Ua(TA,TB) = 0 
%                      (ii): Ub(TA,TB) = 0   [IV].
% However, in order to evaluate the polynomials at the given value pair
% (TA,TB) we need the corresponding coeffient vectors Caa(TB) and Cbb(TA)
% in their quadratic profit-related version.

% Get the vectors of polynomial coefficients.
Caa = CaaPolyQuadratic(C,G,LA,LB,PA,PB,TB,alpha,beta,gamma);
Cbb = CbbPolyQuadratic(C,G,LA,LB,PA,PB,TA,alpha,beta,gamma);
% Evaluate polynomials at the given value pair (TA,TB).
F(1) = polyval(Caa,TA);
F(2) = polyval(Cbb,TB);

% Finally, the nonlinear objective function to be minimized will be given
% as:
%                                    2         2     
%                      F(T) = (F1 - 0)+ (F2 - 0) 
%
% where each Fi with i in {1,2} corresponds to a particular equation of
% the original system.

F = sum(F.^2);

end

