function [RepresentativeSolutions,DigitsAccuracy,ClusterAssignment,Weight] = ExtractClusterRepresentativeSolutions(Solutions)

% This function identifies the "clusters" of almost identical solutions.
% Each "cluster" of "almost identical" solutions is formed by grouping
% together solution istances whose absolute difference is smaller than
% 1.0e-4. In general, a threshold value of the form 1.0e-D indicates that 
% solutions within a give cluster share a sequence of D decimal digits.

% The input argument corresponds to the matrix of the obtained solution
% vectors "Solutions" for which ExitFlag = 1.

% [RepresentativeSolutions] is the output argument that corresponds to the
% identified unique solutions. Each representative solution may be
% interpreted as the center of each cluster but it is actually composed by
% the maximal sequence of identical digits amongts the instances that
% pertain to a certain cluster. RepresentativeSolutions will be a 
% [unique_solutions_number x solution_dimensionality] matrix where
% unique_solutions_number will be the final number of solutions after 
% identifying the different combinations of cluster centers per solution
% dimensionality.

% [DigitsAccuracy] is originally [unique_solutions_number x solution_dimensionality] 
% matrix where each value provides a measure of accuracy of the identified  
% solutions. Therefore, DigitsAccuracy should have a value which is greater  
% or equal to 4 since two solutions are to be considered as identical when  
% at least 4 of their leading digits are identical. Since, we need a single
% accuracy characterization for each unique solution, the DigitsAccuracy
% value will be finally computed by taking the min over all
% solution_dimensions per solution. Therefore, DigitsAccuracy will be a
% [unique_solutions_number x 1] vector. 

% [ClusterAssignment] is a [solutions_number x solution_dimensionality]  
% matrix of cluster indices such that ClusterAssignment(i,j) stores the  
% information concerning the cluster membership of the j-th constituent of  
% the i-th solution. In general, one should extract the unique cluster 
% assignments from matrix ClusterAssignment in order determine the complete 
% set of different solutions. In fact, the unique rows of
% ClusterAssignment matrix determine the number of different solutions
% which will be stored by the variable "unique_solutions_number".

% [Weight] is a [unique_solutions_number x 1] vector indicating the number   
% of instances of "almost identical" solutions that share a maximal sequence 
% of same digits. This value may be interpreted as the "weight" of each 
% representative solution within set of different representative solutions. 
% Ideally, only one representative solution should be obtained. If, on the 
% contrary, more than one representative solutions is retained, after 
% checking the second order optimality conditions, the Weight vector w 
% ill be utilized in order to select the representative solution with the 
% higher weight value.

% Identification of the true number of "clusters" (or representative
% solutions) will be conducted within an iterative procedure where the number
% of clusters will be incrementally changed. The initial number of clusters
% will be 2. Each time the pairwise absolute distances of the associated 
% cluster centers will be larger than the threshold value 1.0e-4 the 
% correponding varible clusters_num will be increased. 

% [ClustersNumber] is an internal [1 x solution_dimension] vector that
% stores the number of clusters that have been identified for each solution
% dimension.

% [ClusterRepresentatives] is an internal [1 x solution_dimension] cell array  
% container which will be storing the representative solutions (or "clusters") 
% identified per solution dimension.

% [MaximumCommonDigits] is an internal [1 x solution_dimension] cell array
% container which will be storing the length of the maximal sequence of
% identical digits for each "cluster" of representative solutions per solution 
% dimension.

% [UniqueClusterAssignment] is an internal [unique_solutions_number x solution_dimensionality]
% matrix which stores the cluster assignment that corresponds to each
% unique representative solution.

% Set the maximum absolute distance.
MaxDist = 1.0e-4; % Actually, the Manhattan distance.

% Get the number of the obtained solutions and the dimensionality of each solution.
[solutions_number,solution_dimensionality] = size(Solutions);

% Initialize the ClusterAssignment matrix.
ClusterAssignment = zeros(solutions_number,solution_dimensionality);

% Initialize the ClustersNumber vector.
ClustersNumber = zeros(1,solution_dimensionality);

% Initialize the ClusterRepresentatives cell array.
ClusterRepresentatives = cell(1,solution_dimensionality);

% Initialize the MaximumCommonDigits cell array.
MaximumCommonDigits  = cell(1,solution_dimensionality);

% Loop through the various dimensions of the solutions matrix.
for dim_index = 1:1:solution_dimensionality
    % Retrieve all solutions for the current solution dimension.
    solution_vector = Solutions(:,dim_index);
    % Initialize a cell array container for the cluster indices per
    % cluster.
    ClusterIndices = cell(1,0);
    % Initialize cell arrays ClusterIndices and ClusterCenters.
    for clusters_num = 1:1:2
        % Get the cluster centers for the current number of clusters.
        [cluster_indices,cluster_centers] = kmeans(solution_vector,clusters_num);
        % Rearrage cluster indices so that they are in accordance with the
        % ascending sorting of the cluster centers.
        [~,ascending_order] = sort(cluster_centers);
        cluster_indices = ascending_order(cluster_indices);
        % Populate cell array ClusterIndices.
        ClusterIndices{clusters_num} = cluster_indices;
        % Get the pairwise manhattan distances between the obtained cluster centers.
        Dc = mandist(cluster_centers,cluster_centers');
        % Get the indices to the diagonal elements of matrix Dc.
        Idiag = [1:clusters_num+1:clusters_num*clusters_num];
        % Set the diagonal elements of Dc equal to 1 so that they do not
        % disturbed the searching process.
        Dc(Idiag) = 1;
    end;
    while(isempty(Dc(Dc<MaxDist)))
        % The number of clusters should be increased.
        clusters_num = clusters_num + 1;
        % Get the cluster centers for the current number of clusters.
        [cluster_indices,cluster_centers] = kmeans(solution_vector,clusters_num);
        % Rearrage cluster indices so that they are in accordance with the
        % ascending sorting of the cluster centers.
        [~,ascending_order] = sort(cluster_centers);
        cluster_indices = ascending_order(cluster_indices);
        % Populate cell array ClusterIndices.
        ClusterIndices{clusters_num} = cluster_indices;
        % Get the pairwise manhattan distances between the obtained cluster centers.
        Dc = mandist(cluster_centers,cluster_centers');
        % Get the indices to the diagonal elements of matrix Dc.
        Idiag = [1:clusters_num+1:clusters_num*clusters_num];
        % Set the diagonal elements of Dc equal to 1 so that they do not
        % disturbed the searching process.
        Dc(Idiag) = 1;
    end;
    % The true number of clusters for the current solution dimension is the
    % one before the last increment.
    clusters_num  = clusters_num - 1;
    % Update the ClustersNumber vector.
    ClustersNumber(dim_index) = clusters_num;
    % Obtain the true cluster indices.
    cluster_indices = ClusterIndices{clusters_num};
    % Update the ClusterAssignment matrix.
    ClusterAssignment(:,dim_index) = cluster_indices;
end;

% Populate the ClusterRepresentatives and MaximumCommonDigits cell arrays.
% Loop through the various solution dimensions.
for dim_index = 1:1:solution_dimensionality
    % Get the number of clusters that were identified for the current
    % solution dimension.
    clusters_num = ClustersNumber(dim_index);
    % Initialize the current element of the ClusterRepresentatives and
    % MaximumCommonDigits cell arrays.
    ClusterRepresentatives{dim_index} = zeros(clusters_num,1);
    MaximumCommonDigits{dim_index} = zeros(clusters_num,1);
    % Loop through the various clusters.
    for cluster_index = 1:1:clusters_num
        % Get all the solution indices for the current solution dimension
        % that pertain to the current cluster.
        current_cluster_indices = find(ClusterAssignment(:,dim_index)==cluster_index);
        % Generate a random pair of indices contained within the current
        % set of cluster indices.
        pair_indices = randi([1 length(current_cluster_indices)],1,2);
        % Get corresponding solutions.
        sol1 = Solutions(current_cluster_indices(pair_indices(1)),dim_index);
        sol2 = Solutions(current_cluster_indices(pair_indices(2)),dim_index);
        % Compute the length of the maximal sequence of identical digits.
        % Convert each numeric expression into a string with 30 decimal
        % points of maximum accuracy. (Most probably there will not be so 
        % many digits)
        csol1 = num2str(sol1,30);
        csol2 = num2str(sol2,30);
        % Get the length of each character sequence.
        l1 = length(csol1);
        l2 = length(csol2);
        % Find the minimum length.
        l_min = min(l1,l2);
        % Adjust the length of each character sequence according to the
        % minimum length.
        csol1 = csol1(1:l_min);
        csol2 = csol2(1:l_min);
        % Compute the vector of element-wise absolute differeces so that
        % same characters will result to the zero value.
        absolute_difference = abs(csol1-csol2);
        % Convert the resulting difference vector into a string so that
        % regular expression routines can be employed.
        absolute_character_difference = num2str(absolute_difference')';
        absolute_character_difference = absolute_character_difference(end,:);
        % Find all sequences of consecutive zeros in the resulting
        % character array.
        zero_tokens = regexp(absolute_character_difference,'0+','match');
        % The required length is the length of the first token minus 2
        % since the first zero and decimal point indicator (0.xxxx).
        max_common_digits = length(zero_tokens{1}) - 2;
        % Compute the value of the representative solution.
        representative_solution = fix(sol1 * (10^max_common_digits)) / (10^max_common_digits + 1);
        % Store representative solution and corresponding legth of the
        % maximal identical digits sequence.
        ClusterRepresentatives{dim_index}(cluster_index) = representative_solution;
        MaximumCommonDigits{dim_index}(cluster_index) = max_common_digits;
    end;
end;

% Extract the unique cluster assignment quadraplets which will finally give
% the total number of different solutions.
UniqueClusterAssignment = unique(ClusterAssignment,'rows');

% Get the number of unique solutions.
unique_solutions_number = size(UniqueClusterAssignment,1);

% Initialize RepresentativeSolutions ans DigitsAccuracy matrices.
RepresentativeSolutions = zeros(unique_solutions_number,solution_dimensionality);
DigitsAccuracy = zeros(unique_solutions_number,solution_dimensionality);

% Compute the weight vector.
% In order to find the absolute frequency of each different solution
% quadraplet we need to transform the indices matrices ClusterAssignment
% and UniqueClusterAssignment into corresponding IDs. Each index quadraplet
% (Q1,Q2,Q3,Q4) will be trasformed to the corresponding decimal number with
% digits D = Q1Q2Q3Q4 and value V(D) = 1000*Q1 + 100*Q2 + 10*Q3 + Q4.

% Set the value vector.
value_vector = 10.^[solution_dimensionality-1:-1:0]';

% Get the cluster assigment matrix IDs.
ClusterAssignmentIDs = ClusterAssignment * value_vector;
% Get the unique cluster assignment matrix IDs.
UniqueClusterAssignmentIDs = UniqueClusterAssignment * value_vector;
% Check whether the number of unique cluster assignment ids is one.
if(length(UniqueClusterAssignmentIDs)==1)
    Weight = [1];
else
    % Get the weight vector.
    Weight = hist(ClusterAssignmentIDs,UniqueClusterAssignmentIDs);
end;

% Populate the RepresentativeSolutions and DigitsAccuracy matrices.
% Loop through the various solition dimensions.
for dim_index = 1:1:solution_dimensionality
    RepresentativeSolutions(:,dim_index) = ClusterRepresentatives{dim_index}(UniqueClusterAssignment(:,dim_index));
    DigitsAccuracy(:,dim_index) = MaximumCommonDigits{dim_index}(UniqueClusterAssignment(:,dim_index));
end;

% Get the final version of the DigitsAccuracy matrix.
DigitsAccuracy = min(DigitsAccuracy,[],2);

end