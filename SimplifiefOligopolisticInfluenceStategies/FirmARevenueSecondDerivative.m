function [DDRa] = FirmARevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the second derivative of the profit for Firm A for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (DDa) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
% Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% DDRa will be composed by the polynomial coefficients and corresponding
% monomial terms of Ha(TA,TB) and Oa(TA,TB) such that:
%             DDRa = Ha / Oa
%             Ha ---> {CHa,THa} and Oa ---> {COa,TOa} 
% where CXa are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {H,O}.

% Set the vectors of multivariate polynomial coefficients CHa and COa.
CHa = [ 4*LA*LB^3*alpha*beta - 4*LA*LB^3*alpha^2 + 4*LA*LB^3*alpha*gamma - ...
        4*LA*LB^3*beta*gamma + 4*C*LA*LB^3*alpha - 4*C*LA*LB^3*beta, 4*LA*LB^4*alpha*gamma ...
        - 4*LA*LB^4*alpha^2 + 4*LA*LB^4*PA*alpha^2 + 4*C*LA*LB^4*alpha - 4*C*LA*LB^4*PA*alpha ...
        - 4*C*LA*LB^4*PB*beta + 4*LA*LB^4*PB*alpha*beta - 4*LA*LB^4*PA*alpha*gamma - 4*LA*LB^4*PB*beta*gamma, ...
        2*LA^2*LB^2*alpha^2 + 6*LA^2*LB^2*beta^2 + 4*C*LA^2*LB^2*alpha - 4*C*LA^2*LB^2*beta ...
        - 8*LA^2*LB^2*alpha*beta + 4*LA^2*LB^2*alpha*gamma - 4*LA^2*LB^2*beta*gamma, ...
        4*LA^2*LB^3*alpha^2 + 8*C*LA^2*LB^3*alpha - 4*C*LA^2*LB^3*beta - 8*LA^2*LB^3*alpha*beta ...
        + 8*LA^2*LB^3*alpha*gamma - 4*LA^2*LB^3*beta*gamma - 8*LA^2*LB^3*PA*alpha^2 ...
        + 12*LA^2*LB^3*PB*beta^2 - 4*C*LA^2*LB^3*PA*alpha - 4*C*LA^2*LB^3*PB*beta ...
        + 12*LA^2*LB^3*PA*alpha*beta - 8*LA^2*LB^3*PB*alpha*beta - ...
        4*LA^2*LB^3*PA*alpha*gamma - 4*LA^2*LB^3*PB*beta*gamma, ...
        2*LA^2*LB^4*alpha^2 + 6*LA^2*LB^4*PA^2*alpha^2 + 6*LA^2*LB^4*PB^2*beta^2 ...
        + 4*C*LA^2*LB^4*alpha + 4*LA^2*LB^4*alpha*gamma - 8*LA^2*LB^4*PA*alpha^2 ...
        - 4*C*LA^2*LB^4*PA*alpha - 4*C*LA^2*LB^4*PB*beta - 8*LA^2*LB^4*PB*alpha*beta ...
        - 4*LA^2*LB^4*PA*alpha*gamma - 4*LA^2*LB^4*PB*beta*gamma + 12*LA^2*LB^4*PA*PB*alpha*beta];

COa = [ LB^4, 4*LA*LB^3, 4*LA*LB^4, 6*LA^2*LB^2, 12*LA^2*LB^3, 6*LA^2*LB^4, ...
        4*LA^3*LB, 12*LA^3*LB^2, 12*LA^3*LB^3, 4*LA^3*LB^4, LA^4, 4*LA^4*LB, ...
        6*LA^4*LB^2, 4*LA^4*LB^3, LA^4*LB^4];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of DDRa for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
  THa = [ TA*TB, TA, TB^2, TB, 1];
  TOa = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^4, TB^3, TB^2, TB, 1];
  Ha = CHa * THa';
  Oa = COa * TOa';
  DDRa = Ha ./ Oa;
% Compute the value of DDRa for the case where TA and TB are matrices of the
% underlying meshgrid.  
else
    THa = { TA.*TB, TA, TB.^2, TB, 1};
    TOa = {TA.^4, (TA.^3).*TB, TA.^3, (TA.^2).*(TB.^2), (TA.^2).*(TB), TA.^2, (TA).*(TB.^3), (TA).*(TB.^2), TA.*TB, TA, TB.^4, TB.^3, TB.^2, TB, 1};
    Ha = zeros(ra,ca);
    Oa = zeros(rb,cb);
    for t = 1:1:length(CHa)
        Ha = Ha + CHa(t)*THa{t};
    end
    for t = 1:1:length(COa)
        Oa = Oa + COa(t)*TOa{t};
    end
    DDRa = Ha ./ Oa;
end

end