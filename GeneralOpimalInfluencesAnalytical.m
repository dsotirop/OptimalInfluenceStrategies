function [T1_opt,T2_opt,Fval,Flag] = GeneralOpimalInfluencesAnalytical(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma)

% This function provides the optimal influences T1_opt, T2_opt for a given
% set of the general optimal influence model parameters that are described  
% by the input variables. Mind tha this function relies exclusively on the
% analytical formulas of the optimal solutions for a given set of input
% arguments and does not utilize any numerical optimizers in order to
% obtained the aforementioned optimal solutions.

% The analytical solutions depend on the roots and sign of the following 
% second degree polynomial:
% G(T1) = (A1-B1)*T^2 - 2*(A0-B0)*T + (A1*B0-A0*B1)

% Output argument Flag denotes whether the solution found is an internal or 
% external solution and the particular type of the exteral solution.
% Flag = 0 : indicates that the solution found was an internal one.
% Flag = 1 : indicates an external minimum solution.
% Flag = 2 : indicates an external maximum solution.

% Compute the coefficients of the polynomial as a function of the model
% parameters.
A2 = -1;
A1 = (Theta1-Theta2) + P2*Lambda2 - P1*Lambda1 + 1;
A0 = Theta2 + P1*(Lambda1 + Theta1*(Lambda1 + Lambda2)) + P2*Theta2*(Lambda1 + Lambda2);
B2 = -1;
B1 = (Lambda2-Lambda1) + (Theta1-Theta2) + 1;
B0 = (Lambda1+Lambda2)*(Theta1+Theta2) + Lambda1 + Theta2;
% Compute the associated discriminant coefficients.
Ga = A1 - B1
Gb = 2*(A0-B0);
Gc = A1*B0 - A0*B1;
DG = Gb^2 - 4 * Ga * Gc;
flag = 0;
flag1 = 0;
flag2 = 0;

% Check whether G(T1) is indeed a second degree polynomial.
if(abs(A1-B1)>eps)
    %warning('G(T1) is a second degree polynomial');
    % Check whether G(T1) has a first order term.
    if(abs(A0-B0)>eps)
        % Check whether G(T)==0 has two roots.
        if(DG>0)
            T1_opt_1 = (-Gb + sqrt(DG))/(2*Ga);
            T1_opt_2 = (-Gb - sqrt(DG))/(2*Ga);
            if(and((T1_opt_1>=Theta1),(T1_opt_1<=1-Theta2)))
                T1_opt = T1_opt_1;
                flag1 = 1;
                Flag = 0;
            end;
            if(and((T1_opt_2>=Theta1),(T1_opt_2<=1-Theta2)))
                T1_opt = T1_opt_2;
                flag2 = 1;
                Flag = 0;
            end;
            if(and(flag1,flag2))
                error('Two valid solutions found');
            end;
            if(and(~flag1,~flag2))
                warning('Both second degree internal solutions are invalid');
            end;
        end;
        % Check whether G(T)==0 has exactly one solution.
        if(DG==0)
            T1_opt = (-Gb)/(2*Ga);
            if(and((T1_opt>=Theta1),(T1_opt<=1-Theta2)))
                flag = 1;
                Flag = 0;
            end;
            if(~flag)
                warning('Single second degree internal solution is invalid');
            end;
        end;
        % Check wheter G(T)==0 has no solutions.
        if(DG<0)
            if(Ga>0)
                T1_opt = 1 - Theta2;
                flag = 1;
                Flag = 2;
            end;
            if(Ga<0)
                T1_opt = Theta1;
                flag = 1;
                Flag = 1;
            end;
        end;
    else
        if(DG>0)
            T1_opt_1 = sqrt(-4*Ga*Gc)/(2*Ga);
            T1_opt_2 = -sqrt(-4*Ga*Gc)/(2*Ga);
            if(and((T1_opt_1>=Theta1),(T1_opt_1<=1-Theta2)))
                T1_opt = T1_opt_1;
                flag1 = 1;
                Flag = 0;
            end;
            if(and((T1_opt_2>=Theta1),(T1_opt_2<=1-Theta2)))
                T1_opt = T1_opt_2;
                flag2 = 1;
                Flag = 0;
            end;
            if(and(flag1,flag2))
                error('Two valid solutions found');
            end;
            if(and(~flag1,~flag2))
                warning('Both second degree internal solutions are invalid');
            end;
        end;
        if(DG==0)
            T1_opt = 0;
            if(and((T1_opt>=Theta1),(T1_opt<=1-Theta2)))
                flag = 1;
                Flag = 0;
            end;
            if(~flag)
                warning('Single second degree internal solution is invalid');
            end;
        end;
        if(DG<0)
            if(Ga>0)
                T1_opt = 1 - Theta2;
                flag = 1;
                Flag = 2;
            end;
            if(Ga<0)
                T1_opt = Theta1;
                flag = 1;
                Flag = 1;
            end;
        end;
    end;
else
    % Check whether G(T1) has a first order term.
    %warning('G(T1) is a first degree polynomial');
    if(abs(A0-B0)>eps)
        T1_opt = -Gc/Gb;
        if(and((T1_opt>=Theta1),(T1_opt<=1-Theta2)))
            flag = 1;
            Flag = 0;
        end;
        if(~flag)
            warning('Single second degree internal solution is invalid');
        end;
    else
        error('The objective function remains constant!')
    end;
end;
if(and((flag==0),and(flag1==0,flag2==0)))
    warning('Entered final if');
    Gmin = Ga*Theta1^2 + Gb*Theta1 + Gc;
    Gmax = Ga*(1-Theta2)^2 + Gb*(1-Theta2) + Gc;
    if(Gmin<0)
        T1_opt = Theta1;
        Flag = 1;
    end;
    if(Gmax>0)
        T1_opt = 1-Theta2;
        Flag = 2;
    end;
    if(and((Gmin<0),(Gmax>0)))
        error('Both outer solutions are valid!')
    end;
end;

if(or(flag,or(flag1,flag2)))
    Gmin = Ga*Theta1^2 + Gb*Theta1 + Gc;
    Gmax = Ga*(1-Theta2)^2 + Gb*(1-Theta2) + Gc;
    if(Gmin<0)
        T1_opt = Theta1;
        Flag = 1;
        warning('Valid solution was found and G(Tmin)<0!');
    end;
    if(Gmax>0)
        T1_opt = 1-Theta2;
        Flag = 2;
        warning('Valid solution was found and G(Tmax)>0!');
    end;
    if(and((Gmin<0),(Gmax>0)))
        error('Both outer solutions are valid!')
    end;
end;
    
T2_opt = 1 - T1_opt;
T_opt = [T1_opt,T2_opt];
Fval = GeneralObjectiveFunction(T_opt,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma);

end

