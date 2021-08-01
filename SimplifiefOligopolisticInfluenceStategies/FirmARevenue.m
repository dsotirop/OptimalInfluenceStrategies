function [Ra] = FirmARevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the revenue for Firm A of the Simplified
% Oligopolistic Optimal Influence Model. The revenue value (Ra) is evaluated
% either on a single pair of (TA,TB) values or on a meshgrid of all possible 
% (TA,TB) pairs which is constructed by corresponding vector of the form 
% tA = [ta_min:dt:ta_max] and tB = [tb_min:dt:tb_max], given the external 
% optimization parameters. Mind that the correct construction of the meshgrid
% implies that the lengths of the initial vectors tA and tB are equal.

% Ra will be composed by the polynomial coefficients and corresponding
% monomial terms of La(TA,TB) and Ma(TA,TB) such that:
%             Ra = La / Ma
%             La ---> {CLa,TLa} and Ma ---> {CMa,TMa} 
% where CXa are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {L,M}.

% Set the vectors of multivariate polynomial coefficients CLa and CMa.
CLa = [ (LB*gamma - LB*alpha + C*LB)^2, 2*(LB*gamma - LB*alpha + C*LB)*(LA*gamma - LA*beta + C*LA), ...
        2*(LB*gamma - LB*alpha + C*LB)*(C*LA*LB + LA*LB*gamma - LA*LB*PA*alpha - LA*LB*PB*beta), ...
        (LA*gamma - LA*beta + C*LA)^2, 2*(LA*gamma - LA*beta + C*LA)*(C*LA*LB + LA*LB*gamma - LA*LB*PA*alpha - LA*LB*PB*beta), ...
        (C*LA*LB + LA*LB*gamma - LA*LB*PA*alpha - LA*LB*PB*beta)^2];

CMa = [ LB^2, 2*LA*LB, 2*LA*LB^2, LA^2, 2*LA^2*LB, LA^2*LB^2];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of Ra for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
    TLa = [ TA^2, TA*TB, TA, TB^2, TB, 1];
    TMa = [ TA^2, TA*TB, TA, TB^2, TB, 1];
    La = CLa * TLa';
    Ma = CMa * TMa';
    Ra = La ./ Ma;
% Compute the value of Ra for the case where TA and TB are matrices of the
% underlying grid.
else
    TLa = { TA.^2, TA.*TB, TA, TB.^2, TB, 1};
    TMa = { TA.^2, TA.*TB, TA, TB.^2, TB, 1};
    La = zeros(ra,ca);
    Ma = zeros(rb,cb);
    for t = 1:1:length(CLa)
        La = La + CLa(t)*TLa{t};
    end
    for t = 1:1:length(CMa)
        Ma = Ma + CMa(t)*TMa{t};
    end
    Ra = La ./ Ma;
end
        

end