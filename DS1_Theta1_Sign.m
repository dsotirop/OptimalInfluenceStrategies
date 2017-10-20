
LAMBDA_RANGE = [0:0.1:0.5];
THETA1_RANGE = [0:0.1:1];
THETA2_RANGE = [0:0.1:1];

Params = [];

for Lambda  = LAMBDA_RANGE
    for Theta1 = THETA1_RANGE
        for Theta2 = THETA2_RANGE
            if(Theta1+Theta2<1)
                DS1 = (2*Lambda*(4*Lambda - 2*Theta1 + 6*Theta2 + 16*Lambda*Theta2 - 2*Theta1*Theta2 - 3*Theta1^2 + 5*Theta2^2 + 1))/(4*Lambda + 2*Theta1 + 2*Theta2 + 8*Lambda*Theta1 + 8*Lambda*Theta2 - 2*Theta1*Theta2 + Theta1^2 + Theta2^2 + 1)^2
                params = [DS1,Lambda,Theta1,Theta2];
                Params = [Params;params];
            end;
        end;
    end;
end;