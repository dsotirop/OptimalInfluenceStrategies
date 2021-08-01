function [RepresentativeSolutions,DigitsAccuracies,RepresentativeFvals] = ExtractClusterRepresentativeSolutions(Solutions,Fvals,MinimumDigitsAccuracy,FirmsNumber)

% This function extracts representative solution tuples for each class of  
% different solutions that are valid within the given set of optimal direct
% influence levels. The given solution tuples are by definition valid since
% they have undergone the filtering process. That is, all available 
% solutions satisfy the First Order and Second Order requirements as well
% as the positivity conditions related to the fundamental quantities of the 
% model. Notice that the optimal Lagrange multiplier solutions will be 
% discarded.
% -------------------------------------------------------------------------
% IMPORTANT NOTE:
% -------------------------------------------------------------------------
% The given set of valid solution points is considered to be the initial
% cluster of solutions. Within this cluster, the maximum pairwise distance
% per solution dimension will be computed. If the maximum of the maximum 
% distances is larger that the maximum distance allowed, then this cluster
% of solutions will be slpit into two additional clusters by utilizing the
% k-means algorithm. Mind that the maximum distance allowed depends upon
% the MinimumDigitsAccuracy parameter which is provided as an input
% argument to the function.
% 
%           (-MinimumDigitsAccuracy + 1)
% Dmax = 10
% 
% The previous condition, is actually the general
% spliting condition of every cluster of solutions that emerges during this
% process. In fact, the indices corresponding to the solution clusters  
% which are generated during this process are assummed to be organized in a 
% binary tree structure which is stored in the cell array ClusterIndexTree.
% The aforementioned cell array will be indexed as:
% SolutionClustersTree{LevelIndex,ClusterIndex} where
% -------------------------------------------------------------------------
% 
%                                        (LevelIndex - 1)
% 1 <= ClusterIndex <= N(LevelIndex) = 2             
% 
% -------------------------------------------------------------------------
% Each time the splitting condition is true for a cluster of solutions
% stored at position SolutionClustersTree{LevelIndex,ClusterIndex} two new
% solution clusters are added to the tree structure stored at positions
% SolutionClustersTree{LevelIndex+1,2*ClusterIndex-1} and
% SolutionClustersTree{LevelIndex+1,2*ClusterIndex}. The final cluster
% representative solutions will be extracted from the leaf level of the
% tree structure if for each cluster of solutions stored in that level the
% splitting condition is found to be false. Mind that a complete row of
% N(ClusterIndex) elements is added to the structure if at least one cluster 
% of solutions at the previous level is found to satisfy the splitting
% condition. Thus, when processing the solution clusters stored at each
% level of the tree structure one has to check whether the current element
% is an empty vector or not.
% -------------------------------------------------------------------------
% Mind that indices corresponding to solution clusters that lie below the 
% root level of the binary tree structure should always be converted to  
% point within the original cluster of solutions stored at:
% SolutionClustersTree{1,1}. 
% -------------------------------------------------------------------------
% Assumming that the variable LevelIndex stores the leaf-level of the
% binary tree structure, then N(LevelIndex) solution clusters must be
% checked at the worst case scenario in order to determine whether the
% iterative procedure for determining the representative solutions must
% continue. Letting Xo the initial boolean value indicating the status of
% the splitting condition before it is evaluated for the first time 
% Xo = 1 and Xj the status of the splitting condition after its evaluation
% at the j-th cluster (1 <= j <= N(LevelIndex)), then the final splitting
% flag Y may be derived by the following equation:
%                                    j=n
%                          /\   \  /
%                 Y = Xo  /  \   \/        Xj
%                                    j=1
% -------------------------------------------------------------------------

% Initialize the output variables as empty matrices.
RepresentativeSolutions = [];
DigitsAccuracies = [];
% Initialize an auxiliary matrix that will store the representative Fvals 
% parameter. In case of multiple solutions, these profit values will be 
% used in order to extract the potentialy unique solution points.
RepresentativeFvals = []; 
% Extract the optimal direct influence level for each firm.
Topt = Solutions(:,1:FirmsNumber);
% Get the number of solutions pertaining to the initial cluster of solution
% points.
SolutionsNumber = size(Solutions,1);
% Initialize the cell array that stores the binary tree structure of the
% solutions clusters.
SolutionClustersTree{1,1} = 1:SolutionsNumber;
% Initialize the index storing the leaf-level of the binary tree.
LevelIndex = 1;
% Initilize the boolean variable indicating whether an extra level has been
% added to the binary tree structure of the solution clusters.
NewLevelFlag = true;

% Start the iterative procedure for determining the representative solutions
% for all different solution tuples present within the initial cluster.
while(NewLevelFlag)
    % Reset the value of the boolean flag indicating the need for adding an
    % extra level at the binary tree structure of solution clusters.
    NewLevelFlag = false;
    % Loop through the indices of the solution clusters that are stored at
    % the leaf-level of the binary stree structure.
    for ClusterIndex = 1:size(SolutionClustersTree,2)
        % Get the indices of solutions pertaining to the current cluster of
        % solutions.
        CurrentNodeIndices = SolutionClustersTree{LevelIndex,ClusterIndex};
        % Check whether the current element of the binary tree structure is 
        % not empty.
        if(~isempty(CurrentNodeIndices))
            % Get the current cluster of solutions.
            CurrentSolutions = Topt(CurrentNodeIndices,:);
            % Get the current cluster of firms' profits.
            CurrentFvals = Fvals(CurrentNodeIndices);
            % Initialize the vector storing the maximum pairwise distance
            % per solution dimension.
            DSOLmax = zeros(1,FirmsNumber);
            % Compute the maximum pairwise distance per solution dimension.
            for firm_index = 1:FirmsNumber
                DSOL = dist(CurrentSolutions(:,firm_index),CurrentSolutions(:,firm_index)');
                DSOL = triu(DSOL);
                DSOL = reshape(DSOL,1,numel(DSOL));
                DSOLmax(firm_index) = max(DSOL);
            end
            % Compute the maximum pairwise distance for both solution
            % dimensions.
            max_dist = max(DSOLmax);
            % If the previously determined splitting condition is false,  
            % then extract a representative solution from the current 
            % cluster.
            if(max_dist <= 10^(-MinimumDigitsAccuracy+1))
                % Update the value of the new tree level boolean flag.
                NewLevelFlag = or(NewLevelFlag,false);
                % Initialize vector containers for storing the representative
                % solution and the corresponding digits accuracy parameter 
                % for the current solution cluster.
                RepresentativeSolution = zeros(1,FirmsNumber);
                DigitsAccuracy = zeros(1,FirmsNumber);
                % Initialize vector container for storing the
                % representative Fvals parameter.
                % the current solution cluster.
                RepresentativeFvals = [RepresentativeFvals,mean(CurrentFvals)];
                % Loop through the various solution dimensions.
                for firm_index = 1:FirmsNumber
                    % Get the representative solution value for the current 
                    % solution dimension.
                    RepresentativeSolution(firm_index) = mean(CurrentSolutions(:,firm_index));
                    % Convert the maximum distance value into a string.
                    % This operation assumes that the maximum distance value is represented
                    % in scientific notation.
                    max_dist_str = num2str(DSOLmax(firm_index));
                    % Get the parts before and after the 'e-xx' segment of the string. 
                    str_parts = strsplit(max_dist_str,'e-');
                    % Check whether the string representation of the max_dist variable
                    % contains indeed two parts.
                    if(length(str_parts)==2)
                        % The string part on the right segment corresponds to the digits
                        % accuracy.
                        DigitsAccuracy(firm_index) = str2double(str_parts{2});
                    else
                        if(max_dist < eps)
                            DigitsAccuracy(firm_index) = Inf;
                        end
                    end
                end
                % Add the representative solution extracted from the
                % current solution cluster to the RepresentativeSolutions
                % matrix.
                RepresentativeSolutions = [RepresentativeSolutions;...
                                           RepresentativeSolution];
                % Add the corresponding digits accuracy parameters to the 
                % DigitsAccuracies matrix.
                DigitsAccuracies = [DigitsAccuracies;DigitsAccuracy];
            % If the previously determined splitting condition is true, 
            % then two new clusters of solutions will be added at the binary
            % tree structure extending its height by one level.
            else
                % Update the value of the new tree level boolean flag.
                NewLevelFlag = or(NewLevelFlag,true);
                % Increase the value of the LevelIndex variable.
                LevelIndex = LevelIndex + 1;
                % Split the current cluster of solutions into two
                % sub-clusters by performing k-means clustering.
                [ClusterIndices,~] = kmeans(CurrentSolutions,2,'Distance',...
                                     'sqeuclidean','Replicates',10,...
                                     'MaxIter',1000);
                % Get the logical array indicating the positions of the
                % solutions that will form the left and right children of 
                % the current node.
                LeftChildIndices = (ClusterIndices==1);
                RightChildIndices = (ClusterIndices==2);
                % Add the left child to the binary tree structure of
                % solution clusters.
                SolutionClustersTree{LevelIndex,2*ClusterIndex-1} = ...
                CurrentNodeIndices(LeftChildIndices);
                % Add the right child to the binary tree structure of
                % solution clusters.
                SolutionClustersTree{LevelIndex,2*ClusterIndex} = ...
                CurrentNodeIndices(RightChildIndices);
            end
        end
    end
end

% Extract unique representative solution or throw a multiple solutions
% error.
[MinFvals,MinIndices] = min(RepresentativeFvals);
% Get the unique min indices.
UniqueMinIndices = unique(MinIndices);
% Check whether a unique representative solution has been found.
if(length(UniqueMinIndices)==1)
    RepresentativeSolutions = RepresentativeSolutions(UniqueMinIndices,:);
    DigitsAccuracies = DigitsAccuracies(UniqueMinIndices,:);
    RepresentativeFvals = MinFvals;
else
% Throw an error indicating the existence of multiple solutions.
    error('Existence of multiple solutions.');
end