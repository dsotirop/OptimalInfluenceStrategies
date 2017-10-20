function Caa = CaaPolyQuadratic(C,G,LA,LB,PA,PB,TB,alpha,beta,gamma)

% This function returns the vector of monomial coefficients for the Caa
% polynomial which corresponds to the reaction function of firm A.

% Mind that this polynomial corresponds to the numerator of the rational
% expression for the first derivative of the profit fa for the first firm
% with respect to TA where fa is assigned a quadratic cost term of the form
% TA^2. The corresponding reaction set is formed by taking into
% consideration only the numerator of the first derivative with respect to
% TA since the associated denumerator is strictly non-zero.

% Define the polynomial coefficient vectors Caa.
Caa = [ -2*G*LB^3, ...
        - 6*G*LA*LB^3 - 6*G*LA*TB*LB^2,...
        - 6*G*LA^2*LB^3 - 12*G*LA^2*LB^2*TB - 6*G*LA^2*LB*TB^2,...
        2*LA*LB^3*alpha^2 - 2*G*LA^3*TB^3 - 2*G*LA^3*LB^3 ...
        - 2*LA*LB^3*alpha*gamma - 6*G*LA^3*LB*TB^2 - 6*G*LA^3*LB^2*TB ...
        - 2*LA*LB^3*PA*alpha^2 + 2*LA*LB^2*TB*alpha^2 - 2*C*LA*LB^3*alpha ...
        + 2*C*LA*LB^3*PA*alpha + 2*C*LA*LB^3*PB*beta - 2*C*LA*LB^2*TB*alpha ...
        + 2*C*LA*LB^2*TB*beta - 2*LA*LB^3*PB*alpha*beta ...
        - 2*LA*LB^2*TB*alpha*beta + 2*LA*LB^3*PA*alpha*gamma ...
        + 2*LA*LB^3*PB*beta*gamma - 2*LA*LB^2*TB*alpha*gamma ...
        + 2*LA*LB^2*TB*beta*gamma, 2*LA^2*LB^3*PA*alpha^2 ...
        - 2*LA^2*LB^3*PB^2*beta^2 - 2*C*LA^2*LB^3*alpha ...
        - 2*LA^2*LB^3*alpha*gamma - 2*LA^2*LB^3*PA^2*alpha^2 ...
        - 2*LA^2*LB*TB^2*beta^2 + 2*LA^2*LB^2*PA*TB*alpha^2 ...
        - 4*LA^2*LB^2*PB*TB*beta^2 + 2*C*LA^2*LB^3*PA*alpha ...
        + 2*C*LA^2*LB^3*PB*beta - 2*C*LA^2*LB*TB^2*alpha ...
        - 4*C*LA^2*LB^2*TB*alpha + 2*C*LA^2*LB*TB^2*beta ...
        + 2*C*LA^2*LB^2*TB*beta + 2*LA^2*LB^3*PB*alpha*beta ...
        + 2*LA^2*LB*TB^2*alpha*beta + 2*LA^2*LB^2*TB*alpha*beta ...
        + 2*LA^2*LB^3*PA*alpha*gamma + 2*LA^2*LB^3*PB*beta*gamma ...
        - 2*LA^2*LB*TB^2*alpha*gamma - 4*LA^2*LB^2*TB*alpha*gamma ...
        + 2*LA^2*LB*TB^2*beta*gamma + 2*LA^2*LB^2*TB*beta*gamma ...
        + 2*C*LA^2*LB^2*PB*TB*beta - 4*LA^2*LB^3*PA*PB*alpha*beta ...
        - 4*LA^2*LB^2*PA*TB*alpha*beta + 2*LA^2*LB^2*PB*TB*alpha*beta ...
        + 2*LA^2*LB^2*PA*TB*alpha*gamma + 2*LA^2*LB^2*PB*TB*beta*gamma ...
        + 2*C*LA^2*LB^2*PA*TB*alpha];


end
