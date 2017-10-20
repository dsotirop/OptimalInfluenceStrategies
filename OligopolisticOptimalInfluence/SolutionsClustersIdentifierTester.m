% This script file provides the main functionality for the algorithmic
% procedure that identifies the "clusters" of similar solutions. Each
% cluster should contain a predefined number of solutions whose values will
% be identical up to a predefined number of decimal points. That is, this
% script generates artificial solutions that form clusters of identical
% values up to randomly chosen decimal point.

% Mind that all solutions lie within the [0,1] interval. Therefore,
% differences can only occur at the decimal points of each number. We 
% assume that the dimensionality of each solution vector is 4.  

clc
clear all

% The number of decimal points that will be identical will vary from
% digit_min to digit_max.
digit_min = 6;
digit_max = 7;

% Set the dimensionality of the solution vectors.
Dimensionality = 4;

% Set the number of clusters to be generated for each solution dimension.
ClustersNum = 3;

% Set the total number of solutions to be generated. 
N = 500;

% Set the cluster indices for each solution per corresponding dimension.
ClusterIndices = randi([1 ClustersNum],N,Dimensionality);

% Generate the main solutions for each cluster per corresponding dimension.
MainSolutions = rand(ClustersNum,Dimensionality);

% Generate the numbers of identical digits per cluster for each solution
% dimension.
Digits = randi([digit_min digit_max],ClustersNum,Dimensionality);

% Initialize the actual solutions matrix.
Solutions = zeros(N,Dimensionality);

% Loop through the various solutions dimensions.
for solution_dimension = 1:1:Dimensionality
    % Loop through the various clusters.
    for cluster_index = 1:1:ClustersNum
        % Find the indices pertaining to the current cluster_index.
        current_cluster_indices = find(ClusterIndices(:,solution_dimension)==cluster_index);
        Solutions(current_cluster_indices,solution_dimension) = MainSolutions(cluster_index,solution_dimension) + rand(1,length(current_cluster_indices)) * 10^(-Digits(cluster_index,solution_dimension)-1);
    end;
end;

% Extract the unique quadraplets of cluster indices.
UniqueClusterIndices = unique(ClusterIndices,'rows');

% Test the matlab routine for the extraction of "clusters" of
% representative solutions.
[RepresentativeSolutions,DigitsAccuracy,ClusterAssignment,Weight] = ExtractClusterRepresentativeSolutions(Solutions)