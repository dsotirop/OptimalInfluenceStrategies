function [DDRb] = FirmBRevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the second derivative of the profit for Firm B for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (DDa) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
% Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% DDRb will be composed by the polynomial coefficients and corresponding
% monomial terms of Hb(TA,TB) and Ob(TA,TB) such that:
%             DDRb = Hb / Ob
%             Hb ---> {CHb,THb} and OB ---> {COb,TOb} 
% where CXb are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {H,O}.

% Set the vectors of multivariate polynomial coefficients CHb and COb.
CHb = [ 2*LA^2*LB^2*alpha^2 + 6*LA^2*LB^2*beta^2 + 4*C*LA^2*LB^2*alpha - ...
        4*C*LA^2*LB^2*beta - 8*LA^2*LB^2*alpha*beta + 4*LA^2*LB^2*alpha*gamma - ...
        4*LA^2*LB^2*beta*gamma, 4*LA^3*LB*alpha*beta - 4*LA^3*LB*alpha^2 + ...
        4*LA^3*LB*alpha*gamma - 4*LA^3*LB*beta*gamma + 4*C*LA^3*LB*alpha - ...
        4*C*LA^3*LB*beta, 4*LA^3*LB^2*alpha^2 + 8*C*LA^3*LB^2*alpha - ...
        4*C*LA^3*LB^2*beta - 8*LA^3*LB^2*alpha*beta + 8*LA^3*LB^2*alpha*gamma - ...
        4*LA^3*LB^2*beta*gamma - 8*LA^3*LB^2*PB*alpha^2 + 12*LA^3*LB^2*PA*beta^2 - ...
        4*C*LA^3*LB^2*PB*alpha - 4*C*LA^3*LB^2*PA*beta - 8*LA^3*LB^2*PA*alpha*beta + ...
        12*LA^3*LB^2*PB*alpha*beta - 4*LA^3*LB^2*PB*alpha*gamma - 4*LA^3*LB^2*PA*beta*gamma, ...
        4*LA^4*LB*alpha*gamma - 4*LA^4*LB*alpha^2 + 4*LA^4*LB*PB*alpha^2 + 4*C*LA^4*LB*alpha - ...
        4*C*LA^4*LB*PB*alpha - 4*C*LA^4*LB*PA*beta + 4*LA^4*LB*PA*alpha*beta - ...
        4*LA^4*LB*PB*alpha*gamma - 4*LA^4*LB*PA*beta*gamma, 2*LA^4*LB^2*alpha^2 + ...
        6*LA^4*LB^2*PB^2*alpha^2 + 6*LA^4*LB^2*PA^2*beta^2 + 4*C*LA^4*LB^2*alpha + ...
        4*LA^4*LB^2*alpha*gamma - 8*LA^4*LB^2*PB*alpha^2 - 4*C*LA^4*LB^2*PB*alpha - ...
        4*C*LA^4*LB^2*PA*beta - 8*LA^4*LB^2*PA*alpha*beta - 4*LA^4*LB^2*PB*alpha*gamma - ...
        4*LA^4*LB^2*PA*beta*gamma + 12*LA^4*LB^2*PA*PB*alpha*beta];

COb = [ LB^4, 4*LA*LB^3, 4*LA*LB^4, 6*LA^2*LB^2, 12*LA^2*LB^3, ...
        6*LA^2*LB^4, 4*LA^3*LB, 12*LA^3*LB^2, 12*LA^3*LB^3, 4*LA^3*LB^4, ...
        LA^4, 4*LA^4*LB, 6*LA^4*LB^2, 4*LA^4*LB^3, LA^4*LB^4];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);    
% Compute the value of DDRb for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
    THb = [ TA^2, TA*TB, TA, TB, 1];
    TOb = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^4, TB^3, TB^2, TB, 1];
    Hb = CHb * THb';
    Ob = COb * TOb';
    DDRb = Hb ./ Ob;
else
% Compute the value of DDRb for the case where TA and TB are matrices of the
% underlying meshgrid.    
    THb = {TA.^2, TA.*TB, TA, TB, 1};
    TOb = {TA.^4, (TA.^3).*(TB), TA.^3, (TA.^2).*(TB.^2), (TA.^2).*(TB), TA.^2, (TA).*(TB.^3), (TA).*(TB.^2), TA.*TB, TA, TB.^4, TB.^3, TB.^2, TB, 1};
    Hb = zeros(ra,ca);
    Ob = zeros(rb,cb);
    for t = 1:1:length(CHb)
        Hb = Hb + CHb(t)*THb{t};
    end
    for t = 1:1:length(COb)
        Ob = Ob + COb(t)*TOb{t};
    end
    DDRb = Hb ./ Ob;
end

end