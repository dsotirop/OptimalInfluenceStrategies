function Cbb = CbbPolyQuadratic(C,G,LA,LB,PA,PB,TA,alpha,beta,gamma)

% This function returns the vector of monomial coefficients for the Cbb
% polynomial which corresponds to the reaction function of firm B.

% Mind that this polynomial corresponds to the numerator of the rational
% expression for the first derivative of the profit fb for the first firm
% with respect to TB where fb is assigned a quadratic cost term of the form
% TB^2. The corresponding reaction set is formed by taking into
% consideration only the numerator of the first derivative with respect to
% TB since the associated denumerator is strictly non-zero.

% Define the polynomial coefficient vectors Cbb.
Cbb = [ -2*G*LA^3, ...
        - 6*G*LB*LA^3 - 6*G*LB*TA*LA^2, ...
        - 6*G*LA^3*LB^2 - 12*G*LA^2*LB^2*TA - 6*G*LA*LB^2*TA^2, ...
        2*LA^3*LB*alpha^2 - 2*G*LB^3*TA^3 - 2*G*LA^3*LB^3 ...
        - 2*LA^3*LB*alpha*gamma - 6*G*LA*LB^3*TA^2 - 6*G*LA^2*LB^3*TA ...
        - 2*LA^3*LB*PB*alpha^2 + 2*LA^2*LB*TA*alpha^2 - 2*C*LA^3*LB*alpha ...
        + 2*C*LA^3*LB*PB*alpha + 2*C*LA^3*LB*PA*beta - 2*C*LA^2*LB*TA*alpha ...
        + 2*C*LA^2*LB*TA*beta - 2*LA^3*LB*PA*alpha*beta ...
        - 2*LA^2*LB*TA*alpha*beta + 2*LA^3*LB*PB*alpha*gamma ...
        + 2*LA^3*LB*PA*beta*gamma - 2*LA^2*LB*TA*alpha*gamma ...
        + 2*LA^2*LB*TA*beta*gamma, 2*LA^3*LB^2*PB*alpha^2 ...
        - 2*LA^3*LB^2*PA^2*beta^2 - 2*C*LA^3*LB^2*alpha ...
        - 2*LA^3*LB^2*alpha*gamma - 2*LA^3*LB^2*PB^2*alpha^2 ...
        - 2*LA*LB^2*TA^2*beta^2 + 2*LA^2*LB^2*PB*TA*alpha^2 ...
        - 4*LA^2*LB^2*PA*TA*beta^2 + 2*C*LA^3*LB^2*PB*alpha ...
        + 2*C*LA^3*LB^2*PA*beta - 2*C*LA*LB^2*TA^2*alpha ...
        - 4*C*LA^2*LB^2*TA*alpha + 2*C*LA*LB^2*TA^2*beta ...
        + 2*C*LA^2*LB^2*TA*beta + 2*LA^3*LB^2*PA*alpha*beta ...
        + 2*LA*LB^2*TA^2*alpha*beta + 2*LA^2*LB^2*TA*alpha*beta ...
        + 2*LA^3*LB^2*PB*alpha*gamma + 2*LA^3*LB^2*PA*beta*gamma ...
        - 2*LA*LB^2*TA^2*alpha*gamma - 4*LA^2*LB^2*TA*alpha*gamma ...
        + 2*LA*LB^2*TA^2*beta*gamma + 2*LA^2*LB^2*TA*beta*gamma ...
        + 2*C*LA^2*LB^2*PA*TA*beta - 4*LA^3*LB^2*PA*PB*alpha*beta ...
        + 2*LA^2*LB^2*PA*TA*alpha*beta - 4*LA^2*LB^2*PB*TA*alpha*beta ...
        + 2*LA^2*LB^2*PB*TA*alpha*gamma + 2*LA^2*LB^2*PA*TA*beta*gamma ...
        + 2*C*LA^2*LB^2*PB*TA*alpha];

end

