function [DRb] = FirmBRevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the first derivative of the revenue for Firm B for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (DRb) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
% Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% DRb will be composed by the polynomial coefficients and corresponding
% monomial terms of Gb(TA,TB) and Sb(TA,TB) such that:
%             DRb = Gb / Sb
%             Gb ---> {CGb,TGb} and Sb ---> {CSb,TSb} 
% where CXb are the multivariate polynomial coefficients with respect to
% both TA and TB and TXb are the corresponding monomial terms, X in {G,S}.

% Set the vectors of multivariate polynomial coefficients CGb and CSb.
CGb = [ -2*(LA*LB*alpha - LA*LB*beta)*(LB*gamma - LB*beta + C*LB), ...
       -2*(LA*LB*alpha - LA*LB*beta)*(LA*gamma - LA*alpha + C*LA), ...
       2*(LA^2*LB*PB*alpha - LA^2*LB*alpha + LA^2*LB*PA*beta)*(LB*gamma - LB*beta + C*LB) ...
       - 2*(LA*LB*alpha - LA*LB*beta)*(C*LA*LB + LA*LB*gamma - LA*LB*PB*alpha - LA*LB*PA*beta), ...
       2*(LA^2*LB*PB*alpha - LA^2*LB*alpha + LA^2*LB*PA*beta)*(LA*gamma - LA*alpha + C*LA), ...
       2*(LA^2*LB*PB*alpha - LA^2*LB*alpha + LA^2*LB*PA*beta)*(C*LA*LB + LA*LB*gamma - ...
       LA*LB*PB*alpha - LA*LB*PA*beta)];
CSb = [ LB^3, 3*LA*LB^2, 3*LA*LB^3, 3*LA^2*LB, 6*LA^2*LB^2, 3*LA^2*LB^3, ...
       LA^3, 3*LA^3*LB, 3*LA^3*LB^2, LA^3*LB^3];
   
% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);  
% Compute the value of DRb for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
    TGb = [ TA^2, TA*TB, TA, TB, 1];
    TSb = [ TA^3, TA^2*TB, TA^2, TA*TB^2, TA*TB, TA, TB^3, TB^2, TB, 1];
    Gb = CGb * TGb';
    Sb = CSb * TSb';
    DRb = Gb ./ Sb;
% Compute the value of DRb for the case where TA and TB are matrices of the
% underlying meshgrid.
else
    TGb = {TA.^2, TA.*TB, TA, TB, 1};
    TSb = {TA.^3, (TA.^2).*TB, TA.^2, (TA).*(TB.^2), TA.*TB, TA, TB.^3, TB.^2, TB, 1};
    Gb = zeros(ra,ca);
    Sb = zeros(rb,cb);
    for t = 1:1:length(CGb)
        Gb = Gb + CGb(t)*TGb{t};
    end
    for t = 1:1:length(CSb)
        Sb = Sb + CSb(t)*TSb{t};
    end
    DRb = Gb ./ Sb;
end

end