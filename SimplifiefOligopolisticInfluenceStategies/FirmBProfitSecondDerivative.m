function [DDb] = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the first derivative of the profit for Firm B for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (DDb) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
%  Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% DDb will be composed by the polynomial coefficients and corresponding
% monomial terms of Wb(TA,TB) and Zb(TA,TB) such that:
%             DDb = Wb / Zb
%             Wb ---> {CWb,TWb} and Zb ---> {CZb,TZb} 
% where CXb are the multivariate polynomial coefficients with respect to
% both TA and TB and TXb are the corresponding monomial terms, X in {W,Z}.

% Set the vectors of multivariate polynomial coefficients CWb and CZb.
CWb = [ -2*G*LB^4, -8*G*LA*LB^3, -8*G*LA*LB^4, -12*G*LA^2*LB^2, ...
        -24*G*LA^2*LB^3, 2*LA^2*LB^2*alpha^2 + 6*LA^2*LB^2*beta^2 - ...
        12*G*LA^2*LB^4 + 4*C*LA^2*LB^2*alpha - 4*C*LA^2*LB^2*beta - ...
        8*LA^2*LB^2*alpha*beta + 4*LA^2*LB^2*alpha*gamma - ...
        4*LA^2*LB^2*beta*gamma, -8*G*LA^3*LB, -24*G*LA^3*LB^2, ...
        4*LA^3*LB*alpha*beta - 24*G*LA^3*LB^3 - 4*LA^3*LB*alpha^2 + ...
        4*LA^3*LB*alpha*gamma - 4*LA^3*LB*beta*gamma + 4*C*LA^3*LB*alpha - ...
        4*C*LA^3*LB*beta, 4*LA^3*LB^2*alpha^2 - 8*G*LA^3*LB^4 + ...
        8*C*LA^3*LB^2*alpha - 4*C*LA^3*LB^2*beta - 8*LA^3*LB^2*alpha*beta + ...
        8*LA^3*LB^2*alpha*gamma - 4*LA^3*LB^2*beta*gamma - ...
        8*LA^3*LB^2*PB*alpha^2 + 12*LA^3*LB^2*PA*beta^2 - ...
        4*C*LA^3*LB^2*PB*alpha - 4*C*LA^3*LB^2*PA*beta - ...
        8*LA^3*LB^2*PA*alpha*beta + 12*LA^3*LB^2*PB*alpha*beta - ...
        4*LA^3*LB^2*PB*alpha*gamma - 4*LA^3*LB^2*PA*beta*gamma, ...
        -2*G*LA^4, -8*G*LA^4*LB, -12*G*LA^4*LB^2, 4*LA^4*LB*alpha*gamma - ...
        8*G*LA^4*LB^3 - 4*LA^4*LB*alpha^2 + 4*LA^4*LB*PB*alpha^2 + ...
        4*C*LA^4*LB*alpha - 4*C*LA^4*LB*PB*alpha - 4*C*LA^4*LB*PA*beta + ...
        4*LA^4*LB*PA*alpha*beta - 4*LA^4*LB*PB*alpha*gamma - ...
        4*LA^4*LB*PA*beta*gamma, 2*LA^4*LB^2*alpha^2 - 2*G*LA^4*LB^4 + ...
        6*LA^4*LB^2*PB^2*alpha^2 + 6*LA^4*LB^2*PA^2*beta^2 + ...
        4*C*LA^4*LB^2*alpha + 4*LA^4*LB^2*alpha*gamma - ...
        8*LA^4*LB^2*PB*alpha^2 - 4*C*LA^4*LB^2*PB*alpha - ...
        4*C*LA^4*LB^2*PA*beta - 8*LA^4*LB^2*PA*alpha*beta - ...
        4*LA^4*LB^2*PB*alpha*gamma - 4*LA^4*LB^2*PA*beta*gamma + ...
        12*LA^4*LB^2*PA*PB*alpha*beta];

CZb = [ LB^4, 4*LA*LB^3, 4*LA*LB^4, 6*LA^2*LB^2, 12*LA^2*LB^3, 6*LA^2*LB^4, ...
        4*LA^3*LB, 12*LA^3*LB^2, 12*LA^3*LB^3, 4*LA^3*LB^4, LA^4, 4*LA^4*LB, ...
        6*LA^4*LB^2, 4*LA^4*LB^3, LA^4*LB^4];    
    
% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of Da for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
   TWb = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^4, TB^3, TB^2, TB, 1];
   TZb = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^4, TB^3, TB^2, TB, 1];
   Wb = CWb * TWb';
   Zb = CZb * TZb';
   DDb = Wb ./ Zb;    
% Compute the value of Da for the case where TA and TB are matrices of the
% underlying meshgrid.
else
   TWb = {TA.^4, (TA.^3).*TB, TA.^3, (TA.^2).*(TB.^2), (TA.^2).*TB, TA.^2, TA.*(TB.^3), TA.*(TB.^2), TA.*TB, TA, TB.^4, TB.^3, TB.^2, TB, 1};
   TZb = {TA.^4, (TA.^3).*TB, TA.^3, (TA.^2).*(TB.^2), (TA.^2).*TB, TA.^2, TA.*(TB.^3), TA.*(TB.^2), TA.*TB, TA, TB.^4, TB.^3, TB.^2, TB, 1};
   Wb = zeros(ra,ca);
   Zb = zeros(rb,cb);
   for t = 1:1:length(CWb)
      Wb = Wb + CWb(t)*TWb{t};      
   end
   for t = 1:1:length(CZb)
      Zb = Zb + CZb(t)*TZb{t};
   end
   DDb = Wb ./ Zb;
end
    
end

