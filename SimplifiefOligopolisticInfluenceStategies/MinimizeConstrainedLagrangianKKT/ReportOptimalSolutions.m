function ReportOptimalSolutions(RSOL,DigitsAccuracy)

% This function reports the optimal solutions along with associated
% DigitsAccuracy and Weight parameters. DigitsAccuracy represents the 
% number of decimal places that were found identical within the set of 
% obtained solutions for each representative solution. The weight parameter
% quantifies the frequency of each representative solution within the set
% of obtained solutions.

% Check whether the set of representative solutions is empty.
if(isempty(RSOL))
    fprintf('No solutions were found meeting the optimization constraints!\n');
else
    % Get the number of representative solutions.
    solutions_number = size(RSOL,1);
    % Inform user in case more than one valid solutions were found.
    if(solutions_number>1)
        fprintf('More than one valid solutions were found!\n');
    end
    % Report the obtained solutions.
    for solution_index = 1:1:solutions_number
        TAopt = RSOL(solution_index,1);
        TBopt = RSOL(solution_index,2);
        digits_accuracy = DigitsAccuracy;
        % weight = Weight(solution_index);
        fprintf('Solution: %d TAopt = %f TBopt = %f\n',solution_index,TAopt,TBopt);
        fprintf('Digits Accuracy for TAopt: %d\n',digits_accuracy(1));
        fprintf('Digits Accuracy for TBopt: %d\n',digits_accuracy(2));
    end

end

