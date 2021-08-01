% This script file provides fundamental visualization functionality for the
% experimentation scenarios stored within the following files:
%
% (i):   'Varying_K.mat'
% (ii):  'Varying_LA_BETA_NEG.mat'
% (iii): 'Varying_PA_BETA_NEG.mat'
% (iv):  'Varying_M.mat'

% Each one from the above .mat files stores an Experiments array of
% Experiment structures.

% Clear command window and workspace.
clc
clear
% Add matlab2tikz folder to the matlab path.
addpath(genpath('matlab2tikz'));
% Set the MinimumDigitsAccuracy parameter used throughout the
% experimentation.
MinimumDigitsAccuracy = 7;
% Set folder containing the experimentation .mat files.
ExperimentFolder = 'experiments';
% Set the filename prefix.
PREFIX = 'Varying_M';
% Set the filename to loaded.
FileName = sprintf('%s.mat',PREFIX);
% Compose the full file name for the current experiment.
CombinedFileName = fullfile(ExperimentFolder,FileName);
% Load the experimentation file.
load(CombinedFileName);
% Determine the number of experiments stored within that particular .mat 
% file. Mind that the number of experiments is always odd.
N = length(Experiments);

% Mind that each experimentation triplet will becomposed by the base line 
% experiment stored at the base_idx position which is always the first 
% Experiment structure store in the array and two more experiments that are
% indexed by plus_idx and minus_idx. Experiment structure indexed by 
% (plus_idx) is diversified from the base line experiment by considering an  
% increased value for one or two of the exogenous parameters that are not 
% allowed to vary. Experiment structure indexed by (minus_idx) is
% diversified from the base line experiment by considering a decreased
% value for one or two of the exogenous parameters that are not allowed to
% vary.

% The term parameters will be identifying the exogenous parameters of the
% model, that is, those variables whose values will not be determined
% through optimization. Thus, the set of parameters will be given as:
%
% Parameters = {LA,LB,PA,PA,K,M,C,G}.
%
% Parameter (G) will be treated as a variable controlling the feasibility
% of the optimization problem.

% The term variables will be identifying the endogenous parameters of the
% model, that is, those variables whose values are determined through
% optimization. Thus, the set of variables will be given as:
%
% Variables = {TA,TB,SA,SB,XA,XB,pA,pB,QA,QB,FA,FB,RA,RB,CA,CB}

% Set the string names of the exogenous parameters in Latex notation.
ParamsNames = {'\lambda_A','\lambda_B','P_A','P_B',...
               '\mu','\kappa','c','\gamma'};
% Get the number of exogenous parameters.
ParamsNum = length(ParamsNames);
% Set the string names of the optimal endogenous parameters to be plotted 
% in Latex notation.
VariablesNames = {'T^{\ast}','s^{\ast}','x^{\ast}','p^{\ast}',...
                  'q^{\ast}','\Pi^{\ast}','R^{\ast}','C^{\ast}'};
% Get the number of endogenous parameters.
VariablesNum = length(VariablesNames);
% Set the figure name strings.
FigureNameStrings = {'Optimal Investment Levels',...
                     'Optimal Limiting Influences',...
                     'Optimal Limiting Beliefs','Optimal Prices',...
                     'Optimal Quantities','Optimal Profits',...
                     'Optimal Revenues','Optimal Costs'};
% Set the file name strings used for saving the plotted figures.
FileNameStrings = {'T_opt','S_opt','X_opt','P_opt','Q_opt','F_opt',...
                   'R_opt','C_opt'};
% Initialize the Triplet structure required by the PlotTriplets and 
% PlotBase functions.
Triplet.PREFIX = PREFIX;
Triplet.FileNameStrings = FileNameStrings;
Triplet.MinimumDigitsAccuracy = MinimumDigitsAccuracy;
Triplet.VariablesNum = VariablesNum;
Triplet.VariablesNames = VariablesNames;
Triplet.ParamsNum = ParamsNum;
Triplet.ParamsNames = ParamsNames;
Triplet.FigureNameStrings = FigureNameStrings;

% Get the number of triplets to be visualized.
Kmax = (N-1) / 2;
% Set the index to the base line experiment of each triplet.
base_idx = 1;
% Load the base line experiment structure.
BaseExperiment = Experiments(base_idx).Experiment;
% Set the base experiment valid indices. Mind that the valid indices are
% those for which the corresponding Topts values are not equal to (-1)
% indicating that a solution point was found for the associated game
% insance.
BaseExperimentValidIndices = (BaseExperiment.Topts(:,1) ~= -1);
% Set the values for the endogenous model parameters for the base line 
% experiment.
BaseExperimentVariables = {BaseExperiment.Topts(BaseExperimentValidIndices,:),...
                           BaseExperiment.Sopts(BaseExperimentValidIndices,:),...
                           BaseExperiment.Xopts(BaseExperimentValidIndices,:),...
                           BaseExperiment.Popts(BaseExperimentValidIndices,:),...
                           BaseExperiment.Qopts(BaseExperimentValidIndices,:),...
                           BaseExperiment.Fopts(BaseExperimentValidIndices,1:2),...
                           BaseExperiment.Fopts(BaseExperimentValidIndices,3:4),...
                           BaseExperiment.Fopts(BaseExperimentValidIndices,5:6)};
% Get the varying parameter index.
VaryingParamIndex = BaseExperiment.ParamIndex;
% Get the varying parameter name.
VaryingParamName = ParamsNames{VaryingParamIndex};
% Get the varying parameter values for the base line experiment.
BaseVaryingParam = BaseExperiment.Parameters{VaryingParamIndex};
% Get the valid varying parameter values for the base line experiment.
BaseVaryingParam = BaseVaryingParam(BaseExperimentValidIndices);
% Get the constant parameters indices.
ConstParamsIndices = setdiff(1:ParamsNum,VaryingParamIndex);
% Get the constant parameter values for the base line experiment.
BaseConstParams = cell2mat(BaseExperiment.Parameters(ConstParamsIndices));
% Get the constant parameter names for the base line experiment.
BaseConstParamsNames = ParamsNames(ConstParamsIndices);

% Update the Triplet structure required by the PlotTriplets and PlotBase 
% functions.
Triplet.VaryingParamName = VaryingParamName;
Triplet.BaseVaryingParam = BaseVaryingParam;
Triplet.BaseExperimentVariables = BaseExperimentVariables;
Triplet.BaseConstParamsNames = BaseConstParamsNames;
Triplet.BaseConstParams = BaseConstParams;

% Plot seperately the base line experiment.
PlotBase(Triplet);
% Plot separately the first derivatives for the base line experiment.
PlotBaseDerivatives(Triplet);

% Loop throught the various triplets.
for k = 1:Kmax
    % Set the index for the plus experiment.
    plus_idx = 2*k;
    % Set the index for the minus experiment.
    minus_idx = 2*k+1;
    % Load the plus experiment structure.
    PlusExperiment = Experiments(plus_idx).Experiment;
    % Set the plus experiment valid indices.
    PlusExperimentValidIndices = (PlusExperiment.Topts(:,1) ~= -1);
    % Get the varying parameter values for the plus experiment.
    PlusVaryingParam = PlusExperiment.Parameters{VaryingParamIndex};
    % Get the valid varying parameter values for the plus experiment.
    PlusVaryingParam = PlusVaryingParam(PlusExperimentValidIndices);
    % Load the minus experiment structure.
    MinusExperiment = Experiments(minus_idx).Experiment;
    % Set the minus experiment valid indices.
    MinusExperimentValidIndices = (MinusExperiment.Topts(:,1) ~= -1);
    % Get the varying parameter values for the minus experiment.
    MinusVaryingParam = MinusExperiment.Parameters{VaryingParamIndex};
    % Get the valid varying parameter values for the minus experiment.
    MinusVaryingParam = MinusVaryingParam(MinusExperimentValidIndices);
    % Get the constant parameter values for the plus experiment.
    PlusConstParams = cell2mat(PlusExperiment.Parameters(ConstParamsIndices));
    % Get the constant parameter values for the minus experiment.
    MinusConstParams = cell2mat(MinusExperiment.Parameters(ConstParamsIndices));
    % Get the constant parameter indices that will be appearing in the
    % title string of each graph.
    TitleConstParamsIndices = (BaseConstParams == PlusConstParams);
    % Express the previously acqired indices with respect to the original
    % indexing of the Parameters cell array.
    TitleConstParamsIndices = ConstParamsIndices(TitleConstParamsIndices);
    % Get the names of the constant parameters that will be appearing on
    % the title of each graph.
    TitleConstParamsNames = ParamsNames(TitleConstParamsIndices);
    % Get the values of the constant parameters that will be appearing on
    % the title of each graph.
    TitleConstParams = cell2mat(BaseExperiment.Parameters(TitleConstParamsIndices));
    % Get the constant parameter indices that will be appearing on the
    % legends of each graph.
    LegendConstParamsIndices = (BaseConstParams ~= PlusConstParams);
    % Express the previously acqired indices with respect to the original
    % indexing of the Parameters cell array.
    LegendConstParamsIndices = ConstParamsIndices(LegendConstParamsIndices);
    % Get the names of the constant parameters that will be appearing on the
    % legends of each graph.
    LegendConstParamsNames = ParamsNames(LegendConstParamsIndices);
    % Get the base line experiment constant parameter values that will be 
    % appearing on the legends of each graph.
    BaseLegendConstParams = cell2mat(BaseExperiment.Parameters(LegendConstParamsIndices));
    % Get the plus experiment constant parameter values that will be
    % appearing on the legends of each graph.
    PlusLegendConstParams = cell2mat(PlusExperiment.Parameters(LegendConstParamsIndices));
    % Get the minus experiment constant parameter values that will be
    % appearing on the legends of each graph.
    MinusLegendConstParams = cell2mat(MinusExperiment.Parameters(LegendConstParamsIndices));
    % Get the plus experiment optimal variables.
    PlusExperimentVariables = {...   
    PlusExperiment.Topts(PlusExperimentValidIndices,:),...        
    PlusExperiment.Sopts(PlusExperimentValidIndices,:),...        
    PlusExperiment.Xopts(PlusExperimentValidIndices,:),...         
    PlusExperiment.Popts(PlusExperimentValidIndices,:),...        
    PlusExperiment.Qopts(PlusExperimentValidIndices,:),...         
    PlusExperiment.Fopts(PlusExperimentValidIndices,1:2),... 
    PlusExperiment.Fopts(PlusExperimentValidIndices,3:4),... 
    PlusExperiment.Fopts(PlusExperimentValidIndices,5:6)};  
    % Set minus experiment variables.
    MinusExperimentVariables = {...
    MinusExperiment.Topts(MinusExperimentValidIndices,:),...        
    MinusExperiment.Sopts(MinusExperimentValidIndices,:),...      
    MinusExperiment.Xopts(MinusExperimentValidIndices,:),...        
    MinusExperiment.Popts(MinusExperimentValidIndices,:),...       
    MinusExperiment.Qopts(MinusExperimentValidIndices,:),...        
    MinusExperiment.Fopts(MinusExperimentValidIndices,1:2),... 
    MinusExperiment.Fopts(MinusExperimentValidIndices,3:4),... 
    MinusExperiment.Fopts(MinusExperimentValidIndices,5:6)};
    
    % Update the Triplet structure required by the PlotTriplets function.
    Triplet.TitleConstParamsNames = TitleConstParamsNames;
    Triplet.TitleConstParams = TitleConstParams;
    Triplet.LegendConstParamsNames = LegendConstParamsNames;
    Triplet.BaseLegendConstParams = BaseLegendConstParams;
    Triplet.PlusLegendConstParams = PlusLegendConstParams;
    Triplet.MinusLegendConstParams = MinusLegendConstParams;
    Triplet.PlusExperimentVariables = PlusExperimentVariables;
    Triplet.MinusExperimentVariables = MinusExperimentVariables;
    Triplet.PlusVaryingParam = PlusVaryingParam;
    Triplet.MinusVaryingParam = MinusVaryingParam;
    
    % Plot model variables for the current experimental triplet.
    PlotTriplets(Triplet,k);
end

% Define a function for plotting the triplets of model variables with
% respect to the varying parameter.
function PlotTriplets(Triplet,triplet_index)

% Triplet: is a structure storing all the necessary information for carrying
%          out required plots.
% triplex_index: identifies the current experimentation triplet index.     

% Define the rgb colors to be used for plotting the variables of each
% experiment.
BlackColor =  [0,0,0];
RedColor   =  [(255/255),(0/255),(0/255)];
GreenColor =  [(0/255),(128/255),(0/255)];
BlueColor  =  [(0/255),(0/255),(255/255)];
DeepSkyBlueColor = [(0/255),(191/255),(255/255)];

% Set color the Consumer plots.
color_c = DeepSkyBlueColor;

% Loop through the various model variables that need to be plotted.
% Mind that the optimal limiting influences, optimal revenues and optimal
% costs are not ploted.
for variable_index = [1,2,3,5,6]
    % Create a new plotting figure.
    f = figure('Name',Triplet.FigureNameStrings{variable_index});
    % Set the current variable name.
    VariableName = Triplet.VariablesNames{variable_index};
    % Set the current varying parameter name.
    VaryingParamName = Triplet.VaryingParamName;
    % Set the x-label string for the current figure.
    XLabelString  = VaryingParamName;
    % Set the y-label string for the current figure.
    YLabelString = sprintf('%s_A / %s_B',VariableName,VariableName);
    % Set the first line of the title string for the current figure.
    if(variable_index ~= 2)
        FirstLineTitleString = sprintf('%s_A(%s) vs %s_B(%s)',VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    else
       FirstLineTitleString = sprintf('%s_A(%s) vs %s_B(%s) vs %s_C(%s)',VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    end
    % Set the second line of the title string for the current figure.
    L = length(Triplet.TitleConstParams);
    for t = 1:L
        if(t==1)
            SecondLineTitleString = sprintf('%s = %s',...
                                        Triplet.TitleConstParamsNames{t},...
                                        num2str(Triplet.TitleConstParams(t)));
        else
            SecondLineTitleString = sprintf('%s  %s = %s',...
                                        SecondLineTitleString,...
                                        Triplet.TitleConstParamsNames{t},...
                                        num2str(Triplet.TitleConstParams(t)));
        end
    end
    % Set the third , fourth and fifth lines of the title string denoting
    % the parameter diversification considered in the minus, base and plus
    % experiments.
    ThirdLineTitleString  = '(i):';
    FourthLineTitleString = '(ii):';
    FifthLineTitleString  = '(iii):';
    Lo = length(Triplet.LegendConstParamsNames);
    for t = 1:Lo
        if(t==Lo)
            ThirdLineTitleString = sprintf('%s %s = %s',...
                                   ThirdLineTitleString,...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.MinusLegendConstParams(t)));
            FourthLineTitleString = sprintf('%s %s = %s',...
                                    FourthLineTitleString,...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.BaseLegendConstParams(t)));
            FifthLineTitleString = sprintf('%s %s = %s',...
                                   FifthLineTitleString,...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.PlusLegendConstParams(t)));                    
        else
            ThirdLineTitleString = sprintf('%s %s = %s ',...
                                   ThirdLineTitleString,...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.MinusLegendConstParams(t)));
            FourthLineTitleString = sprintf('%s %s = %s ',...
                                    FourthLineTitleString,...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.BaseLegendConstParams(t)));
            FifthLineTitleString = sprintf('%s %s = %s ',...
                                   FifthLineTitleString,...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.PlusLegendConstParams(t)));                               
        end
    end
    % Use one or five lines of the title string.
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString};
    %TitleString = {FirstLineTitleString,SecondLineTitleString,...
                   %ThirdLineTitleString,FourthLineTitleString,...
                   %FifthLineTitleString};
    % Set the title for the current figure.
    title(TitleString,'fontsize',14);
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',14);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',14);
    
    % Plot the optimal variables for the two firms with respect to the
    % varying parameter. Mind that when the optimal limiting influences are
    % plot an additional plot for the limiting influence of the consumer
    % should be incorporated in the graph.
    hold on
    % Set the independent variables of each plot appearing on the x-axis.
    x_minus = Triplet.MinusVaryingParam;
    x_base = Triplet.BaseVaryingParam;
    x_plus = Triplet.PlusVaryingParam;
    % Set the dependent variables of each plot appearing on the y-axis.
    y_minus = Triplet.MinusExperimentVariables{variable_index};
    y_base = Triplet.BaseExperimentVariables{variable_index};
    y_plus = Triplet.PlusExperimentVariables{variable_index};
    % Set the colors for the curves representing the optimal variables for
    % the two firms.
    if(variable_index ~= 2)
        dy_minus = abs(y_minus(:,1)-y_minus(:,2));
        dy_base = abs(y_base(:,1)-y_base(:,2));
        dy_plus = abs(y_plus(:,1)-y_plus(:,2));
    else
        dy_minus = abs(y_minus(:,1)-y_minus(:,3));
        dy_base = abs(y_base(:,1)-y_base(:,3));
        dy_plus = abs(y_plus(:,1)-y_plus(:,3));
    end
    dy_minus = (1/length(dy_minus))*sum(dy_minus);
    dy_base = (1/length(dy_base))*sum(dy_base);
    dy_plus = (1/length(dy_plus))*sum(dy_plus);
    if(dy_minus < 10^(-Triplet.MinimumDigitsAccuracy))
        color_a_minus = BlueColor;
        color_b_minus = BlueColor;
    else
        color_a_minus = RedColor;
        color_b_minus = GreenColor;
    end
    if(dy_base < 10^(-Triplet.MinimumDigitsAccuracy))
        color_a_base = BlueColor;
        color_b_base = BlueColor;
    else
        color_a_base = RedColor;
        color_b_base = GreenColor;
    end
    if(dy_plus < 10^(-Triplet.MinimumDigitsAccuracy))
        color_a_plus = BlueColor;
        color_b_plus = BlueColor;
    else
        color_a_plus = RedColor;
        color_b_plus = GreenColor;
    end
    % Do the actual plotting for the minus experiment for the two firms.
    if(variable_index ~= 2)
        p_minus_a = plot(x_minus,y_minus(:,1),'-.','LineWidth',2.3);
        p_minus_b = plot(x_minus,y_minus(:,2),'-.','LineWidth',2.3);
        set(p_minus_a,'Color',color_a_minus);
        set(p_minus_b,'Color',color_b_minus);
    else
        p_minus_a = plot(x_minus,y_minus(:,1),'-.','LineWidth',2.3);
        p_minus_b = plot(x_minus,y_minus(:,3),'-.','LineWidth',2.3);
        set(p_minus_a,'Color',color_a_minus);
        set(p_minus_b,'Color',color_b_minus);
        p_minus_c = plot(x_minus,y_minus(:,2),'-.','LineWidth',2.3);
        set(p_minus_c,'Color',color_c);
    end
    % Do the actual plotting for the base experiment for the two firms.
    if(variable_index ~= 2)
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
    else
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,3),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
        p_base_c = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_c,'Color',color_c);
    end
    % Do the actual plotting for the plus experiment for the two firms.
    if(variable_index ~= 2)
        p_plus_a = plot(x_plus,y_plus(:,1),':','LineWidth',2.3);
        p_plus_b = plot(x_plus,y_plus(:,2),':','LineWidth',2.3);
        set(p_plus_a,'Color',color_a_plus);
        set(p_plus_b,'Color',color_b_plus);
    else
        p_plus_a = plot(x_plus,y_plus(:,1),':','LineWidth',2.3);
        p_plus_b = plot(x_plus,y_plus(:,3),':','LineWidth',2.3);
        set(p_plus_a,'Color',color_a_plus);
        set(p_plus_b,'Color',color_b_plus);
        p_plus_c = plot(x_plus,y_plus(:,2),':','LineWidth',2.3);
        set(p_plus_c,'Color',color_c);
    end
    % Enable plotting grid.
    grid on
    % Set the axes limits.
    x_min = min([min(x_minus) min(x_base) min(x_plus)]);
    x_max = max([max(x_minus) max(x_base) max(x_plus)]);
    y_min = min([min(min(y_minus)) min(min(y_base)) min(min(y_plus))]);
    y_max = max([max(max(y_minus)) max(min(y_base)) max(max(y_plus))]);
    if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
        y_min = 0;
        y_max = 1;
    end
    axis([x_min x_max y_min y_max]);
    % Set the linewidth and color of the current axes.
    h = gca;
    h.LineWidth = 2.5;
    set(h,'XColor',BlackColor);
    set(h,'YColor',BlackColor);
    hold off
    % Set the size and location of the plotting figure.
    set(f,'Position',[100 100 500 400])
    
    % Save the current plotting figure in eps format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.eps']));
    print(figurefilename,'-depsc','-r2000','-noui');
   
    % Save the current plotting figure in jpg format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.jpg']));
    print(figurefilename,'-djpeg','-r1000','-noui');
    
    % Save the currebt ploting figure in pdf format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.pdf']));
    print(figurefilename,'-dpdf','-r2000','-noui');
    
    % Set the current plotting figure in tex format.
    % texfilename = fullfile('figures',Triplet.PREFIX,...
    % strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.tex']));
    % matlab2tikz(texfilename);
end
end

% Define a function for plotting the base experiment's model variables with
% respect to the varying parameter.
function PlotBase(Triplet)

% Triplet: is a structure storing all the necessary information for carrying
%          out required plots. Mind that when the function PlotBase is called 
%          the Triplet structure may not yet contain entries regarding the
%          plus and minus experiment.

% Define the rgb colors to be used for plotting the variables of each
% experiment.
BlackColor =  [0,0,0];
RedColor   =  [(255/255),(0/255),(0/255)];
GreenColor =  [(0/255),(128/255),(0/255)];
BlueColor  =  [(0/255),(0/255),(255/255)];
OrangeColor = [(255/255),(165/255),(0/255)];

% Set color the Consumer plots.
color_c = OrangeColor;

% Loop through the various model variables that need to be plotted.
for variable_index = 1:Triplet.VariablesNum
    % Create a new plotting figure.
    f = figure('Name',Triplet.FigureNameStrings{variable_index});
    % Set the current variable name.
    VariableName = Triplet.VariablesNames{variable_index};
    % Set the current varying parameter name.
    VaryingParamName = Triplet.VaryingParamName;
    % Set the x-label string for the current figure.
    XLabelString  = VaryingParamName;
    % Set the y-label string for the current figure.
    YLabelString = sprintf('$$\\mathbf{%s_A / %s_B}$$',VariableName,VariableName);
    % Set the first line of the title string for the current figure.
    if(variable_index ~= 2)
        FirstLineTitleString = sprintf('$$\\mathbf{%s_A(%s)~\\textrm{vs}~%s_B(%s)}$$',VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    else
       % Set the y-label string for the current figure.
       YLabelString = sprintf('$$\\mathbf{%s_A / %s_B / %s_C}$$',VariableName,VariableName,VariableName);
       FirstLineTitleString = sprintf('$$\\mathbf{%s_A(%s)~\\textrm{vs}~%s_B(%s)~\\textrm{vs}~%s_C(%s)}$$',...
                                       VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    end
    % Set the second line of the title string for the current figure.
    L = Triplet.ParamsNum - 1;
    for t = 1:L
        if(t==1)
            SecondLineTitleString = sprintf('%s = %s',...
                                        Triplet.BaseConstParamsNames{t},...
                                        num2str(Triplet.BaseConstParams(t)));
        else
            SecondLineTitleString = sprintf('%s  %s = %s',...
                                        SecondLineTitleString,...
                                        Triplet.BaseConstParamsNames{t},...
                                        num2str(Triplet.BaseConstParams(t)));
        end
    end
    % Use one or two lines of the title string.
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString};
    % TitleString = {FirstLineTitleString,SecondLineTitleString};
    
    % Set the title for the current figure.
    title(TitleString,'fontsize',14,'Interpreter','latex');
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',14);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',10,'Interpreter','latex');
    
    % Plot the optimal variables for the two firms with respect to the
    % varying parameter. Mind that when the optimal limiting influences are
    % plot an additional plot for the limiting influence of the consumer
    % should be incorporated in the graph.
    hold on
    % Set the independent variables of each plot appearing on the x-axis.
    x_base = Triplet.BaseVaryingParam;
    % Set the dependent variables of each plot appearing on the y-axis.
    y_base = Triplet.BaseExperimentVariables{variable_index};
    % Set the colors for the curves representing the optimal variables for
    % the two firms.
    if(variable_index ~= 2)
        dy_base = abs(y_base(:,1)-y_base(:,2));
    else
        dy_base = abs(y_base(:,1)-y_base(:,3));
    end
    dy_base = (1/length(dy_base))*sum(dy_base);
    if(dy_base < 10^(-Triplet.MinimumDigitsAccuracy))
        color_a_base = BlueColor;
        color_b_base = BlueColor;
    else
        color_a_base = RedColor;
        color_b_base = GreenColor;
    end
    % Do the actual plotting for the base experiment for the two firms.
    if(variable_index ~= 2)
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
    else
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,3),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
        p_base_c = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_c,'Color',color_c);
    end
    % Enable plotting grid.
    grid on
    % Set the axes limits.
    x_min = min(x_base);
    x_max = max(x_base);
    y_min = min(min(y_base));
    y_max = max(max(y_base));
    if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
        y_min = 0;
        y_max = 1;
    end
    axis([x_min x_max y_min y_max]);
    % Set the linewidth and color of the current axes.
    h = gca;
    h.LineWidth = 2.5;
    set(h,'XColor',BlackColor);
    set(h,'YColor',BlackColor);
    hold off
    % Set the size and location of the plotting figure.
    set(f,'Position',[100 100 500 400])
    
    % Save the current plotting figure in eps format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_0.eps']));
    print(figurefilename,'-depsc','-r2000','-noui');
   
    % Save the current plotting figure in jpg format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_0.jpg']));
    print(figurefilename,'-djpeg','-r1000','-noui');
    
    % Save the current plotting figure in pdf format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_0.pdf']));
    print(figurefilename,'-dpdf','-r2000','-noui');
    
    % Set the current plotting figure in tex format.
    % texfilename = fullfile('figures',Triplet.PREFIX,...
    % strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.tex']));
    % matlab2tikz(texfilename);
end

end

% Define a function for plotting the base experiment's model variables with
% respect to the varying parameter.
function PlotBaseDerivatives(Triplet)

% Triplet: is a structure storing all the necessary information for carrying
%          out required plots. Mind that when the function PlotBase is called 
%          the Triplet structure may not yet contain entries regarding the
%          plus and minus experiment.

% Define the rgb colors to be used for plotting the variables of each
% experiment.
BlackColor =  [0,0,0];
RedColor   =  [(255/255),(0/255),(0/255)];
GreenColor =  [(0/255),(128/255),(0/255)];
BlueColor  =  [(0/255),(0/255),(255/255)];
OrangeColor = [(255/255),(165/255),(0/255)];

% Set color the Consumer plots.
color_c = OrangeColor;

% Loop through the various model variables that need to be plotted.
for variable_index = 1:Triplet.VariablesNum
    % Create a new plotting figure.
    f = figure('Name',Triplet.FigureNameStrings{variable_index});
    % Set the current variable name.
    VariableName = Triplet.VariablesNames{variable_index};
    % Set the current varying parameter name.
    VaryingParamName = Triplet.VaryingParamName;
    % Set the x-label string for the current figure.
    XLabelString  = VaryingParamName;
    % Set the y-label string for the current figure.
    YLabelString = sprintf('$$\\mathbf{\\frac{d%s_A}{d%s} / \\frac{d%s_B}{d%s}}$$',VariableName,VaryingParamName,...
                           VariableName,VaryingParamName);
    % Set the first line of the title string for the current figure.
    if(variable_index ~= 2)
        FirstLineTitleString = sprintf('$$\\mathbf{\\frac{d%s_A}{d%s}~\\textrm{vs}~\\frac{d%s_B}{d%s}}$$',...
                                       VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    else
       YLabelString = sprintf('$$\\mathbf{\\frac{d%s_A}{d%s} / \\frac{d%s_B}{d%s} / \\frac{d%s_C}{d%s}}$$',...
                                       VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
       FirstLineTitleString = sprintf('$$\\mathbf{\\frac{d%s_A}{d%s}~\\textrm{vs}~\\frac{d%s_B}{d%s}~\\textrm{vs}~\\frac{d%s_C}{d%s}}$$',...
                                       VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName,VariableName,...
                                       VaryingParamName);
    end
    % Set the second line of the title string for the current figure.
    L = Triplet.ParamsNum - 1;
    for t = 1:L
        if(t==1)
            SecondLineTitleString = sprintf('%s = %s',...
                                        Triplet.BaseConstParamsNames{t},...
                                        num2str(Triplet.BaseConstParams(t)));
        else
            SecondLineTitleString = sprintf('%s  %s = %s',...
                                        SecondLineTitleString,...
                                        Triplet.BaseConstParamsNames{t},...
                                        num2str(Triplet.BaseConstParams(t)));
        end
    end
    % Use one or two lines of the title string.
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString};
    %TitleString = {FirstLineTitleString,SecondLineTitleString};
    
    % Set the title for the current figure.
    title(TitleString,'fontsize',14,'Interpreter','latex');
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',14);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',10,'Interpreter','latex');

    % Plot the optimal variables for the two firms with respect to the
    % varying parameter. Mind that when the optimal limiting influences are
    % plot an additional plot for the limiting influence of the consumer
    % should be incorporated in the graph.
    hold on
    % Set the independent variables of each plot appearing on the x-axis.
    x_base = Triplet.BaseVaryingParam;
    % Set the dependent variables of each plot appearing on the y-axis.
    y_base = Triplet.BaseExperimentVariables{variable_index};
    % Compute the derivatives of the dependent variables.
    y_base = diff(y_base) ./ diff(x_base)';
    % Adjust the length of the independent variable to reflect the range of
    % values of the derivatives of the dependent variables.
    x_base = x_base(2:end);
    % Set the colors for the curves representing the optimal variables for
    % the two firms.
    if(variable_index ~= 2)
        dy_base = abs(y_base(:,1)-y_base(:,2));
    else
        dy_base = abs(y_base(:,1)-y_base(:,3));
    end
    dy_base = (1/length(dy_base))*sum(dy_base);
    if(dy_base < 10^(-Triplet.MinimumDigitsAccuracy))
        color_a_base = BlueColor;
        color_b_base = BlueColor;
    else
        color_a_base = RedColor;
        color_b_base = GreenColor;
    end
    % Do the actual plotting for the base experiment for the two firms.
    if(variable_index ~= 2)
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
    else
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.3);
        p_base_b = plot(x_base,y_base(:,3),'-','LineWidth',2.3);
        set(p_base_a,'Color',color_a_base);
        set(p_base_b,'Color',color_b_base);
        p_base_c = plot(x_base,y_base(:,2),'-','LineWidth',2.3);
        set(p_base_c,'Color',color_c);
    end
    % Enable plotting grid.
    grid on
    % Set the axes limits.
    x_min = min(x_base);
    x_max = max(x_base);
    y_min = min(min(y_base));
    y_max = max(max(y_base));
    if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
        y_min = 0;
        y_max = 1;
    end
    axis([x_min x_max y_min y_max]);
    % Set the linewidth and color of the current axes.
    h = gca;
    h.LineWidth = 2.5;
    set(h,'XColor',BlackColor);
    set(h,'YColor',BlackColor);
    hold off
    % Set the size and location of the plotting figure.
    set(f,'Position',[100 100 500 400])
    
    % Save the current plotting figure in eps format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_00.eps']));
    print(figurefilename,'-depsc','-r2000','-noui');
   
    % Save the current plotting figure in jpg format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_00.jpg']));
    print(figurefilename,'-djpeg','-r1000','-noui');
    
    % Save the current plotting figure in pdf format.
    figurefilename = fullfile('figures',Triplet.PREFIX,...
    strcat([Triplet.FileNameStrings{variable_index} '_00.pdf']));
    print(figurefilename,'-dpdf','-r2000','-noui');
    
    % Set the current plotting figure in tex format.
    % texfilename = fullfile('figures',Triplet.PREFIX,...
    % strcat([Triplet.FileNameStrings{variable_index} '_' num2str(triplet_index) '.tex']));
    % matlab2tikz(texfilename);
end

end