function [Rb] = FirmBRevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the revenue for Firm B of the Simplified
% Oligopolistic Optimal Influence Model. The revenue value (Rb) is evaluated
% either on a single pair of (TA,TB) values or on a meshgrid of all possible 
% (TA,TB) pairs which is constructed by corresponding vector of the form 
% tA = [ta_min:dt:ta_max] and tB = [tb_min:dt:tb_max], given the external 
% optimization parameters. Mind that the correct construction of the meshgrid
% implies that the lengths of the initial vectors tA and tB are equal.

% Ra will be composed by the polynomial coefficients and corresponding
% monomial terms of Lb(TA,TB) and Mb(TA,TB) such that:
%             Rb = Lb / Mb
%             Lb ---> {CLb,TLb} and Mb ---> {CMb,TMb} 
% where CXb are the multivariate polynomial coefficients with respect to
% both TA and TB and TXb are the corresponding monomial terms, X in {L,M}.

% Set the vectors of multivariate polynomial coefficients CLb and CMb.
CLb = [ (LB*gamma - LB*beta + C*LB)^2, 2*(LA*gamma - LA*alpha + C*LA)*(LB*gamma - LB*beta + C*LB),...
        2*(LB*gamma - LB*beta + C*LB)*(C*LA*LB + LA*LB*gamma - LA*LB*PB*alpha - LA*LB*PA*beta), ...
        (LA*gamma - LA*alpha + C*LA)^2, 2*(LA*gamma - LA*alpha + C*LA)*(C*LA*LB + LA*LB*gamma - LA*LB*PB*alpha - LA*LB*PA*beta), ...
        (C*LA*LB + LA*LB*gamma - LA*LB*PB*alpha - LA*LB*PA*beta)^2];
    
CMb = [ LB^2, 2*LA*LB, 2*LA*LB^2, LA^2, 2*LA^2*LB, LA^2*LB^2];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of Rb for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
    TLb = [ TA^2, TA*TB, TA, TB^2, TB, 1];
    TMb = [ TA^2, TA*TB, TA, TB^2, TB, 1];
    Lb = CLb * TLb';
    Mb = CMb * TMb';
    Rb = Lb ./ Mb;
% Compute the value of Rb for the case where TA and TB are matrices of the
% underlying grid.
else
    TLb = { TA.^2, TA.*TB, TA, TB.^2, TB, 1};
    TMb = { TA.^2, TA.*TB, TA, TB.^2, TB, 1};
    Lb = zeros(ra,ca);
    Mb = zeros(rb,cb);
    for t = 1:1:length(CLb)
        Lb = Lb + CLb(t)*TLb{t};
    end
    for t = 1:1:length(CMb)
        Mb = Mb + CMb(t)*TMb{t};
    end
    Rb = Lb ./ Mb;
end

end