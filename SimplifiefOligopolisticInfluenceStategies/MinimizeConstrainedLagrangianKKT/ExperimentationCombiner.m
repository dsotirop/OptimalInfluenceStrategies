% This script file combines within a single structure the experimentation
% results obtained for a specific exogenous parameter. The different
% experimentation scenarios correspond to the following mat files:
%
% (1): 'Varying_K_X.mat'(where X in {1,...,10})
% (2): 'Varying_LA_BETA_NEG_X.mat' (where X in {1,...,12})
% (3): 'Varying_PA_BETA_NEG_X.mat' (where X in {1,...,12})
% (4): 'Varying_M_X.mat' (where X in {1,...,10})

% Clear command window and workspace.
clc
clear
% Set folder containing the experimentation .mat files.
ExperimentFolder = 'experiments';
% Set the string description of the experimentation scenario.
Description = 'Varying_M';
% Set the minimum and maximum indices for the experimentation datasets 
% matching the previously set description.
min_experiment_idx = 0;
max_experiment_idx = 6;
% Initialize placeholder.
ExperimentsCell = cell(1,max_experiment_idx-max_experiment_idx+1);
% Set the filename of the combined experimentation structure for the
% current description.
CombinedFileName = fullfile(ExperimentFolder,Description);
% Loop through the various experimentation datasets for the current
% description.
for k = min_experiment_idx:max_experiment_idx
    % Set the name of the current experimentation file to be loaded.
    current_filename = sprintf('%s_%d.mat',Description,k);
    % Set the full filename for the current experiment.
    current_full_filename = fullfile(ExperimentFolder,current_filename);
    % Load .mat file.
    load(current_full_filename);
    % Append the current experimentation structure to the array.
    ExperimentsCell{k+1} = Experiment;
    % Clear variable.
    clear Experiment 
end

% Convert cell array of structures to structures array
Experiments = cell2struct(ExperimentsCell,'Experiment',1);
% Save the overall structure.
save(CombinedFileName,'Experiments');
