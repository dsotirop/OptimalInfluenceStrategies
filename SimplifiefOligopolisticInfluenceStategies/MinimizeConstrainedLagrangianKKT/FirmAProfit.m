function [fa] = FirmAProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function evaluates the profit for Firm A of the Simplified
% Oligopolistic Optimal Influence Model. The profit value (fa) is evaluated
% either on a single pair of (TA,TB) values or on a meshgrid of all possible 
% (TA,TB) pairs which is constructed by corresponding vector of the form 
% tA = [ta_min:dt:ta_max] and tB = [tb_min:dt:tb_max], given the external 
% optimization parameters. Mind that the correct construction of the meshgrid
% implies that the lengths of the initial vectors tA and tB are equal.

% fa will be composed by the polynomial coefficients and corresponding
% monomial terms of Pa(TA,TB) and Qa(TA,TB) such that:
%             fa = Pa / Qa
%             Pa ---> {CPa,TPa} and Qa ---> {CQa,TQa} 
% where CXa are the multivariate polynomial coefficients with respect to
% both TA and TB and TXa are the corresponding monomial terms, X in {P,Q}.

% Set the vectors of multivariate polynomial coefficients CPa and CQa.
CPa = [ -G*LB^2, -2*G*LA*LB, -2*G*LA*LB^2, -G*LA^2, -2*G*LA^2*LB, ...
        C^2*LB^2 - 2*C*LB^2*alpha + 2*C*LB^2*gamma - G*LA^2*LB^2 + ...
        LB^2*alpha^2 - 2*LB^2*alpha*gamma + LB^2*gamma^2, ...
        2*C^2*LA*LB + 2*LA*LB*gamma^2 - 2*C*LA*LB*alpha - 2*C*LA*LB*beta + ...
        4*C*LA*LB*gamma + 2*LA*LB*alpha*beta - 2*LA*LB*alpha*gamma - ...
        2*LA*LB*beta*gamma, 2*LA*LB^2*gamma^2 + 2*C^2*LA*LB^2 + ...
        4*C*LA*LB^2*gamma - 2*LA*LB^2*alpha*gamma + 2*LA*LB^2*PA*alpha^2 - ...
        2*C*LA*LB^2*alpha - 2*C*LA*LB^2*PA*alpha - 2*C*LA*LB^2*PB*beta + ...
        2*LA*LB^2*PB*alpha*beta - 2*LA*LB^2*PA*alpha*gamma - ...
        2*LA*LB^2*PB*beta*gamma, ...
        C^2*LA^2 - 2*C*LA^2*beta + 2*C*LA^2*gamma + LA^2*beta^2 - ...
        2*LA^2*beta*gamma + LA^2*gamma^2, ....
        2*LA^2*LB*gamma^2 + 2*C^2*LA^2*LB + 4*C*LA^2*LB*gamma - ...
        2*LA^2*LB*beta*gamma + 2*LA^2*LB*PB*beta^2 - 2*C*LA^2*LB*beta - ...
        2*C*LA^2*LB*PA*alpha - 2*C*LA^2*LB*PB*beta + 2*LA^2*LB*PA*alpha*beta - ...
        2*LA^2*LB*PA*alpha*gamma - 2*LA^2*LB*PB*beta*gamma, ...
        C^2*LA^2*LB^2 - 2*C*LA^2*LB^2*PA*alpha - 2*C*LA^2*LB^2*PB*beta + ...
        2*C*LA^2*LB^2*gamma + LA^2*LB^2*PA^2*alpha^2 + ...
        2*LA^2*LB^2*PA*PB*alpha*beta - 2*LA^2*LB^2*PA*alpha*gamma + ...
        LA^2*LB^2*PB^2*beta^2 - 2*LA^2*LB^2*PB*beta*gamma + LA^2*LB^2*gamma^2];

CQa = [ LB^2, 2*LA*LB, 2*LA*LB^2, LA^2, 2*LA^2*LB, LA^2*LB^2];

% Get the sizes of input variable TA and TB.
[ra,ca] = size(TA);
[rb,cb] = size(TB);
% Compute the value of fa for the case where TA and TB are single-valued
% vectors.
if(ra*ca*rb*cb==1)
   TPa = [ TA^4, TA^3*TB, TA^3, TA^2*TB^2, TA^2*TB, TA^2, TA*TB, TA, TB^2, TB, 1];
   TQa = [ TA^2, TA*TB, TA, TB^2, TB, 1];
   Pa = CPa * TPa';
   Qa = CQa * TQa';
   fa = Pa ./ Qa;    
% Compute the value of fa for the case where TA and TB are matrices of the
% underlying meshgrid.
else
   TPa = { TA.^4, (TA.^3).*TB, TA.^3, (TA.^2).*(TB.^2), (TA.^2).*TB, TA.^2, TA.*TB, TA, TB.^2, TB, 1};
   TQa = { TA.^2, TA.*TB, TA, TB.^2, TB, 1};
   Pa = zeros(ra,ca);
   Qa = zeros(rb,cb);
   for t = 1:1:length(CPa)
      Pa = Pa + CPa(t)*TPa{t};      
   end;
   for t = 1:1:length(CQa)
      Qa = Qa + CQa(t)*TQa{t};
   end;
   fa = Pa ./ Qa;
end;


end
