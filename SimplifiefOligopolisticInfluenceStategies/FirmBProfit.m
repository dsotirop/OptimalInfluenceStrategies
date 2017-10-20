function [fb] = FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the profit for Firm B of the Simplified
% Oligopolistic Optimal Influence Model. The profit value (fb) is evaluated
% either on a single pair of (TA,TB) values or on a meshgrid of all possible 
% (TA,TB) pairs which is constructed by corresponding vector of the form 
% tA = [ta_min:dt:ta_max] and tB = [tb_min:dt:tb_max], given the external 
% optimization parameters. Mind that the correct construction of the meshgrid
% implies that the lengths of the initial vectors tA and tB are equal.

% fa will be composed by the polynomial coefficients and corresponding
% monomial terms of Pb(TA,TB) and Qb(TA,TB) such that:
%             fb = Pb / Qb
%             Pb ---> {CPb,TPb} and Qb ---> {CQb,TQb} 
% where CXb are the multivariate polynomial coefficients with respect to
% both TA and TB and TXb are the corresponding monomial terms, X in {P,Q}.

% Set the vectors of multivariate polynomial coefficients CPb and CQb.
CPb = [ -G*LB^2, ...
        C^2*LB^2 - 2*C*LB^2*beta + 2*C*LB^2*gamma + LB^2*beta^2 - ...
        2*LB^2*beta*gamma + LB^2*gamma^2, -2*G*LA*LB, -2*G*LA*LB^2, ...
        2*C^2*LA*LB + 2*LA*LB*gamma^2 - 2*C*LA*LB*alpha - 2*C*LA*LB*beta + ...
        4*C*LA*LB*gamma + 2*LA*LB*alpha*beta - 2*LA*LB*alpha*gamma - ...
        2*LA*LB*beta*gamma, 2*LA*LB^2*gamma^2 + 2*C^2*LA*LB^2 + ...
        4*C*LA*LB^2*gamma - 2*LA*LB^2*beta*gamma + 2*LA*LB^2*PA*beta^2 - ...
        2*C*LA*LB^2*beta - 2*C*LA*LB^2*PB*alpha - 2*C*LA*LB^2*PA*beta + ...
        2*LA*LB^2*PB*alpha*beta - 2*LA*LB^2*PB*alpha*gamma - ...
        2*LA*LB^2*PA*beta*gamma, -G*LA^2, -2*G*LA^2*LB, ...
        C^2*LA^2 - 2*C*LA^2*alpha + 2*C*LA^2*gamma - G*LA^2*LB^2 + ...
        LA^2*alpha^2 - 2*LA^2*alpha*gamma + LA^2*gamma^2, ...
        2*LA^2*LB*gamma^2 + 2*C^2*LA^2*LB + 4*C*LA^2*LB*gamma - ...
        2*LA^2*LB*alpha*gamma + 2*LA^2*LB*PB*alpha^2 - 2*C*LA^2*LB*alpha ...
        - 2*C*LA^2*LB*PB*alpha - 2*C*LA^2*LB*PA*beta + ...
        2*LA^2*LB*PA*alpha*beta - 2*LA^2*LB*PB*alpha*gamma - ...
        2*LA^2*LB*PA*beta*gamma, C^2*LA^2*LB^2 - 2*C*LA^2*LB^2*PA*beta - ...
        2*C*LA^2*LB^2*PB*alpha + 2*C*LA^2*LB^2*gamma + LA^2*LB^2*PA^2*beta^2 + ...
        2*LA^2*LB^2*PA*PB*alpha*beta - 2*LA^2*LB^2*PA*beta*gamma + ...
        LA^2*LB^2*PB^2*alpha^2 - 2*LA^2*LB^2*PB*alpha*gamma + LA^2*LB^2*gamma^2];

CQb = [ LB^2, 2*LA*LB, 2*LA*LB^2, LA^2, 2*LA^2*LB, LA^2*LB^2];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of fb for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
   TPb = [ TA^2*TB^2, TA^2, TA*TB^3, TA*TB^2, TA*TB, TA, TB^4, TB^3, TB^2, TB, 1];
   TQb = [ TA^2, TA*TB, TA, TB^2, TB, 1];
   Pb = CPb * TPb';
   Qb = CQb * TQb';
   fb = Pb ./ Qb;    
% Compute the value of fb for the case where TA and TB are matrices of the
% underlying meshgrid.
else
   TPb = { (TA.^2).*(TB.^2), TA.^2, TA.*(TB.^3), TA.*(TB.^2), TA.*TB, TA, TB.^4, TB.^3, TB.^2, TB, 1};
   TQb = {TA.^2, TA.*TB, TA, TB.^2, TB, 1};
   Pb = zeros(ra,ca);
   Qb = zeros(rb,cb);
   for t = 1:1:length(CPb)
      Pb = Pb + CPb(t)*TPb{t};      
   end;
   for t = 1:1:length(CQb)
      Qb = Qb + CQb(t)*TQb{t};
   end;
   fb = Pb ./ Qb;
end;

end

