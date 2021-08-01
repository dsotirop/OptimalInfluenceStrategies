function [Da] = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the first derivative of the profit for Firm A for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (Da) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
%  Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% Da will be composed by the polynomial coefficients and corresponding
% monomial terms of Ua(TA,TB) and Va(TA,TB) such that:
%             Da = Ua / Va
%             Ua ---> {CUa,TVa} and Qa ---> {CVa,TVa} 
% where CXa are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {U,V}.

% Set the vectors of multivariate polynomial coefficients CUa and CVa.
CUa = [ -2*G*LB^3, -6*G*LA*LB^2, -6*G*LA*LB^3, -6*G*LA^2*LB, ...
        -12*G*LA^2*LB^2, -6*G*LA^2*LB^3, -2*G*LA^3, -6*G*LA^3*LB, ...
        2*LA*LB^2*alpha^2 - 6*G*LA^3*LB^2 - 2*LA*LB^2*alpha*beta - ...
        2*LA*LB^2*alpha*gamma + 2*LA*LB^2*beta*gamma - 2*C*LA*LB^2*alpha + ...
        2*C*LA*LB^2*beta, 2*LA*LB^3*alpha^2 - 2*G*LA^3*LB^3 - ...
        2*LA*LB^3*alpha*gamma - 2*LA*LB^3*PA*alpha^2 - 2*C*LA*LB^3*alpha + ...
        2*C*LA*LB^3*PA*alpha + 2*C*LA*LB^3*PB*beta - 2*LA*LB^3*PB*alpha*beta + ...
        2*LA*LB^3*PA*alpha*gamma + 2*LA*LB^3*PB*beta*gamma, ...
        2*LA^2*LB*alpha*beta - 2*LA^2*LB*beta^2 - 2*LA^2*LB*alpha*gamma + ...
        2*LA^2*LB*beta*gamma - 2*C*LA^2*LB*alpha + 2*C*LA^2*LB*beta, ...
        2*C*LA^2*LB^2*beta - 4*C*LA^2*LB^2*alpha + 2*LA^2*LB^2*alpha*beta ...
        - 4*LA^2*LB^2*alpha*gamma + 2*LA^2*LB^2*beta*gamma + ...
        2*LA^2*LB^2*PA*alpha^2 - 4*LA^2*LB^2*PB*beta^2 + ...
        2*C*LA^2*LB^2*PA*alpha + 2*C*LA^2*LB^2*PB*beta - ...
        4*LA^2*LB^2*PA*alpha*beta + 2*LA^2*LB^2*PB*alpha*beta + ...
        2*LA^2*LB^2*PA*alpha*gamma + 2*LA^2*LB^2*PB*beta*gamma, ...
        2*LA^2*LB^3*PA*alpha^2 - 2*LA^2*LB^3*PB^2*beta^2 - ...
        2*C*LA^2*LB^3*alpha - 2*LA^2*LB^3*alpha*gamma - ...
        2*LA^2*LB^3*PA^2*alpha^2 + 2*C*LA^2*LB^3*PA*alpha + ...
        2*C*LA^2*LB^3*PB*beta + 2*LA^2*LB^3*PB*alpha*beta + ...
        2*LA^2*LB^3*PA*alpha*gamma + 2*LA^2*LB^3*PB*beta*gamma - ...
        4*LA^2*LB^3*PA*PB*alpha*beta];

CVa = [ LB^3, 3*LA*LB^2, 3*LA*LB^3, 3*LA^2*LB, 6*LA^2*LB^2, ...
        3*LA^2*LB^3, LA^3, 3*LA^3*LB, 3*LA^3*LB^2, LA^3*LB^3];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of Da for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
   TUa = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^2, TB, 1];
   TVa = [ TA^3, TA^2*TB, TA^2, TA*TB^2, TA*TB, TA, TB^3, TB^2, TB, 1];
   Ua = CUa * TUa';
   Va = CVa * TVa';
   Da = Ua ./ Va;    
% Compute the value of Da for the case where TA and TB are matrices of the
% underlying meshgrid.
else
   TUa = { TA.^4, (TA.^3).*TB, TA.^3, (TA.^2).*(TB.^2), (TA.^2).*TB, TA.^2, TA.*(TB.^3), TA.*(TB.^2), TA.*TB, TA, TB.^2, TB, 1};
   TVa = { TA.^3, (TA.^2).*TB, TA.^2, TA.*(TB.^2), TA.*TB, TA, TB.^3, TB.^2, TB, 1};
   Ua = zeros(ra,ca);
   Va = zeros(rb,cb);
   for t = 1:1:length(CUa)
      Ua = Ua + CUa(t)*TUa{t};      
   end
   for t = 1:1:length(CVa)
      Va = Va + CVa(t)*TVa{t};
   end
   Da = Ua ./ Va;
end

end

