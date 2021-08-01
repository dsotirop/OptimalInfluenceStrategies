 function [RepresentativeSolution,DigitsAccuracy] = ExtractRepresentativeSolution(Solutions,MinimumDigitsAccuracy)

% This function extracts a representative solution from a given set of
% solutions under the assumptions that the obtained solutions have
% undergone the filtering process. That is, all available solutions satisfy
% the First Order and Second Order requirements as well as the positivity 
% conditions related to the fundamental quantities of the model. Notice  
% that the optimal Lagrange multiplier solutions have been discarded.
% -------------------------------------------------------------------------
% IMPORTANT NOTE:
% -------------------------------------------------------------------------
% It has to be mentioned that this function is unable to extract multiple
% representative solutions in case such solutions do exist and satisfy all
% the necessary first and second order conditions. However, an error will 
% be raised in case the extracted representative solution violates the
% minimum digits accuracy directive. Such a violation might be indicative
% of the non-uniqueness of the extracted solution.
% -------------------------------------------------------------------------

% Get the dimensionality of the solution vectors.
Dimensionality = size(Solutions,2);
% Initialize a vector storing the maximum pairwise distance per solution
% dimension.
DSOL_max = zeros(1,Dimensionality);
% Initialize a vector storing the digits accuracy per solution dimension.
DigitsAccuracy = zeros(1,Dimensionality);
% Initalize a vector storing the representative solution value per solution
% dimension.
RepresentativeSolution = zeros(1,Dimensionality);

% Loop through the various dimensions of the solution vector.
for dim_index = 1:Dimensionality
    % Get the pairwise distances for all the solution points.
    DSOL = dist(Solutions(:,dim_index),Solutions(:,dim_index)');
    % Get the upper triangular part of the distance matrix and reshape it into
    % a row vector.
    DSOL = triu(DSOL);
    DSOL = reshape(DSOL,1,numel(DSOL));
    % Get the maxixum pairwise distance.
    max_dist = max(DSOL);
    DSOL_max(dim_index) = max_dist;
    % Check whether the maximum pairwise distance exceeds the maximum 
    % allowed distance based on the minimum digits accuracy parameter.
    if(max_dist > 10^(-MinimumDigitsAccuracy+1))
        error('Accuracy Violation: dimension index = %d digits_accuracy = %d\n',dim_index,DigitsAccuracy(dim_index));
    else
        % Get the representative solution value for the current solution
        % dimension.
        RepresentativeSolution(dim_index) = mean(Solutions(:,dim_index));
        % Convert the maximum distance value into a string.
        % This operation assumes that the maximum distance value is represented
        % in scientific notation.
        max_dist_str = num2str(max_dist);
        % Get the parts before and after the 'e-xx' segment of the string. 
        str_parts = strsplit(max_dist_str,'e-');
        % Check whether the string representation of the max_dist variable
        % contains indeed two parts.
        if(length(str_parts)==2)
            % The string part on the right segment corresponds to the digits
            % accuracy.
            DigitsAccuracy(dim_index) = str2double(str_parts{2});
        else
            if(max_dist < eps)
                DigitsAccuracy(dim_index) = Inf;
            end
        end
    end
end

end