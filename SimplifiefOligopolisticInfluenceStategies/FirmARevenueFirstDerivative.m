function [DRa] = FirmARevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the first derivative of the revenue for Firm A for 
% the Simplified Oligopolistic Optimal Influence Model. The derivative value 
% (DRa) is evaluated either on a single pair of (TA,TB) values or on a 
% meshgrid of all possible  (TA,TB) pairs which is constructed by  
% corresponding vector of the form tA = [ta_min:dt:ta_max] and  
% tB = [tb_min:dt:tb_max], given the external optimization parameters. 
% Mind that the correct construction of the meshgrid implies that the 
% lengths of the initial vectors tA and tB are equal.

% DRa will be composed by the polynomial coefficients and corresponding
% monomial terms of Ga(TA,TB) and Sa(TA,TB) such that:
%             DRa = Ga / Sa
%             Ga ---> {CGa,TGa} and Sa ---> {CSa,TSa} 
% where CXa are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {G,S}.

% Set the vectors of multivariate polynomial coefficients CGa and CSa.
CGa = [ -2*(LA*LB*alpha - LA*LB*beta)*(LB*gamma - LB*alpha + C*LB), ...
        2*(LA*LB^2*PA*alpha - LA*LB^2*alpha + LA*LB^2*PB*beta)*(LB*gamma - LB*alpha + C*LB), ...
        -2*(LA*LB*alpha - LA*LB*beta)*(LA*gamma - LA*beta + C*LA), ...
        2*(LA*LB^2*PA*alpha - LA*LB^2*alpha + LA*LB^2*PB*beta)*(LA*gamma - LA*beta + C*LA) - ...
        2*(LA*LB*alpha - LA*LB*beta)*(C*LA*LB + LA*LB*gamma - LA*LB*PA*alpha - LA*LB*PB*beta), ...
        2*(LA*LB^2*PA*alpha - LA*LB^2*alpha + LA*LB^2*PB*beta)*(C*LA*LB + LA*LB*gamma - ...
        LA*LB*PA*alpha - LA*LB*PB*beta)];
CSa = [ LB^3, 3*LA*LB^2, 3*LA*LB^3, 3*LA^2*LB, 6*LA^2*LB^2, 3*LA^2*LB^3, ...
        LA^3, 3*LA^3*LB, 3*LA^3*LB^2, LA^3*LB^3]; 

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);    
% Compute the value of DRa for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
    TGa = [ TA*TB, TA, TB^2, TB, 1];
    TSa = [ TA^3, TA^2*TB, TA^2, TA*TB^2, TA*TB, TA, TB^3, TB^2, TB, 1];
    Ga = CGa * TGa';
    Sa = CSa * TSa';
    DRa = Ga ./ Sa;
else
% Compute the value of DRa for the case where TA and TB are matrices of the
% underlying meshgrid.    
    TGa = { TA.*TB, TA, TB.^2, TB, 1};
    TSa = { TA.^3, (TA.^2).*(TB), TA.^2, (TA).*(TB.^2), TA.*TB, TA, TB.^3, TB.^2, TB, 1};
    Ga = zeros(ra,ca);
    Sa = zeros(rb,cb);
    for t = 1:1:length(CGa)
        Ga = Ga + CGa(t)*TGa{t};
    end
    for t = 1:1:length(CSa)
        Sa = Sa + CSa(t)*TSa{t};
    end
    DRa = Ga ./ Sa;
end


end