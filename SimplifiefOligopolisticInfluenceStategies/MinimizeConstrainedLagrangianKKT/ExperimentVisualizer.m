% This script file provides fundamental visualization functionality for a
% particular experimentation scenario stored within one of the following
% files:
%
%   (i): 'Varying_PA_BETA_NEG_X.mat' with X in {0,...,12}
%  (ii): 'Varying_LA_BETA_NEG_X.mat' with X in {0,...,12}
% (iii): 'Varying_K_X.mat' with X in {0,...,6}
%  (iv): 'Varying_M_X.mat' with X in {0,...,6}

% Each one from the above .mat files stores an Experiment structure storing
% the individual information concerning the endogenous and exogenous model
% variables for the corresponding experiment.

% Clear command window and workspace.
clc
clear
% Set the MinimumDigitsAccuracy parameter used throughout the
% experimentation.
MinimumDigitsAccuracy = 7;
% Set folder containing the experimentation .mat files.
ExperimentFolder = 'experiments';
% Set the filename prefix.
PREFIX = 'Varying_K_6';
% Set the filename to loaded.
FileName = sprintf('%s.mat',PREFIX);
% Compose the full file name for the current experiment.
CombinedFileName = fullfile(ExperimentFolder,FileName);
% Load the experimentation file.
load(CombinedFileName);

% Determine the valid ids for which the optimimation process terminated
% successfully.
idx = find(Experiment.Topts(:,1)~=-1);

% Get the values of the varying parameter (t) for the current experiment.
t = Experiment.Parameters{Experiment.ParamIndex};
% Get the name of the varying parameter.
ParamName = Experiment.ParamsNames{Experiment.ParamIndex};
% Get the valid values for the varying parameter.
t = t(idx)';
% Get the internal parameters LA, LB, mu and kappa.
LA = Experiment.Parameters{1};
if(length(LA)>1)
    LA = LA(idx)';
end
LB = Experiment.Parameters{2};
if(length(LB)>1)
    LB = LB(idx)';
end
mu = Experiment.Parameters{5};
if(length(mu)>1)
    mu = mu(idx)';
end
kappa = Experiment.Parameters{6};
if(length(kappa)>1)
    kappa = kappa(idx)';
end
% Compute the values of the parameters alpha and beta.
alpha = (kappa .* mu - 2) ./ (mu.^2 - 4);
beta = (2*kappa - mu) ./ (mu.^2 - 4);

% Get the optimal investment levels per firm.
TA = Experiment.Topts(idx,1);
TB = Experiment.Topts(idx,2);
% Compute the values of the quantity So = LA*LB + LB*TA + LA*TB.
So = LA .* LB + LB .* TA + LA .* TB;
% Compute an approximation for the first derivatives of the optimization
% variables TA and TB with respect to the varying exogenous parameter t.
dTA = diff(TA);
dTB = diff(TB);
dt = diff(t);
dTAdt = dTA ./ dt;
dTBdt = dTB ./ dt;
% Define a structure of additional input arguments to be used by the
% plotting functions.
ExtraParams.LA = LA;
ExtraParams.LB = LB;
ExtraParams.alpha = alpha;
ExtraParams.beta = beta;
ExtraParams.So = So;
ExtraParams.dTAdt = dTAdt;
ExtraParams.dTBdt = dTBdt;
% Set the current variable index.
VariableIndex = 1;
% Plot curves TA and TB along with their first derivatives.
PlotCurves(VariableIndex,Experiment.ParamIndex,TA,TB,t,MinimumDigitsAccuracy);

% Get the optimal limiting influences per firm.
SA = Experiment.Sopts(idx,1);
SB = Experiment.Sopts(idx,3);
% Compute an approximation for the first derivatives of the optimization
% variables SA and SB with respect to the varying exogenous parameter t.
dSA = diff(SA);
dSB = diff(SB);
dt = diff(t);
dSAdt = dSA ./ dt;
dSBdt = dSB ./ dt;
% Set the current variable index.
VariableIndex = 2;
% Plot curves SA and SB along with their first derivatives.
PlotCurves(VariableIndex,Experiment.ParamIndex,SA,SB,t,MinimumDigitsAccuracy);

% Get the optimal limiting beliefs per firm.
XA = Experiment.Xopts(idx,1);
XB = Experiment.Xopts(idx,2);
% Compute an approximation for the first derivatives of the optimization
% variables XA and XB with respect to the varying exogenous parameter t.
dXA = diff(XA);
dXB = diff(XB);
dt = diff(t);
dXAdt = dXA ./ dt;
dXBdt = dXB ./ dt;
% Set the current variable index.
VariableIndex = 3;
% Plot curves SA and SB along with their first derivatives.
PlotCurves(VariableIndex,Experiment.ParamIndex,XA,XB,t,MinimumDigitsAccuracy);

% Get the optimal quantities per firm.
QA = Experiment.Qopts(idx,1);
QB = Experiment.Qopts(idx,2);
% Compute an approximation for the first derivatives of the optimization
% variables QA and QB with respect to the varying exogenous parameter t.
dQA = diff(QA);
dQB = diff(QB);
dt = diff(t);
dQAdt = dQA ./ dt;
dQBdt = dQB ./ dt;
% Set the current variable index.
VariableIndex = 5;
% Plot curves QA and QB along with their first derivatives.
PlotCurves(VariableIndex,Experiment.ParamIndex,QA,QB,t,MinimumDigitsAccuracy);

% Get the optimal profits per firm.
FA = Experiment.Fopts(idx,1);
FB = Experiment.Fopts(idx,4);
% Set the current variable index.
VariableIndex = 6;
% Plot curves FA and FB along with their first derivatives.
PlotCurves(VariableIndex,Experiment.ParamIndex,FA,FB,t,MinimumDigitsAccuracy);

function [color_a,color_b] = PlotColors(VA,VB,MinimumDigitsAccuracy)

% This function checks whether the optimal values for variables VA and VB
% are the same so as to use the same color for plotting purposes.

% Define the rgb colors to be used for plotting the variables of each
% experiment.
RedColor   =  [(255/255),(0/255),(0/255)];
GreenColor =  [(0/255),(128/255),(0/255)];
BlueColor  =  [(0/255),(0/255),(255/255)];

% Get the absolute difference of values between variables VA and VB.
dV = abs(VA-VB);
% Compute the average absolute difference.
dV = (1/length(dV)) * sum(dV);

% If the average absolute difference is less than the minimum digits
% accuracy, then variables VA and VB are considered to have the same
% values.
if(dV < 10^(-MinimumDigitsAccuracy))
    color_a = BlueColor;
    color_b = BlueColor;
else
    color_a = RedColor;
    color_b = GreenColor;
end

end

function PlotCurves(VariableIndex,ParamIndex,VA,VB,t,MinimumDigitsAccuracy)

% This function performs the actual plotting of curves VA and VB as well as
% the first derivatives of the quantities VA and VB with respect to the
% varying parameter.

% Set the string names of the optimal endogenous parameters to be plotted 
% in Latex notation.
VariablesNames = {'T^{\ast}','s^{\ast}','x^{\ast}','p^{\ast}',...
                  'q^{\ast}','\Pi^{\ast}','R^{\ast}','C^{\ast}'};

% Set the string names of the exogenous parameters in Latex notation.
ParamsNames = {'\lambda_A','\lambda_B','p_A^{(0)}','p_B^{(0)}',...
               '\mu','\kappa','c','\gamma'};

% Set the figure name strings.
FigureNameStrings = {'Optimal Investment Levels',...
                     'Optimal Limiting Influences',...
                     'Optimal Limiting Beliefs',...
                     'Optimal Prices',...
                     'Optimal Quantities',...
                     'Optimal Profits',...
                     'Optimal Revenues','Optimal Costs'};
                 
% Set the derivatives figure name strings.
DerivativeFigureNameStrings = {'Optimal Investment Levels Derivatives',...
                               'Optimal Limiting Influences Derivatives',...
                               'Optimal Limiting Beliefs Derivatives',...
                               'Optimal Prices Derivatives',...
                               'Optimal Quantities Derivatives',...
                               'Optimal Profits Derivatives',...
                               'Optimal Revenues Derivatives',...
                               'Optimal Costs Derivatives'};

% Set the name of the endogenous optimization variable to be plotted.
VariableName = VariablesNames{VariableIndex};
% Set the name of the exogenous model parameter to be plotted.
ParamName = ParamsNames{ParamIndex};
% Set the name of the figure.
FigureName = FigureNameStrings{VariableIndex};
% Set the title of the figure.
FigureTitle = sprintf('%s_A(%s) vs %s_B(%s)',VariableName,ParamName,VariableName,ParamName);
% Set the x-axis label of the figure.
XLabel = ParamName;
% Set the y-axis lable of the figure.
YLabel = sprintf('%s_A / %s_B',VariableName,VariableName);
% Get the color for each curve.
[color_a,color_b] = PlotColors(VA,VB,MinimumDigitsAccuracy);

% Do the actual plotting for the optimal optimization variables of the 
% model with respect to the varying exogenous parameter.
figure('Name',FigureName);
title(FigureTitle,'fontsize',14);
xlabel(XLabel,'fontsize',14);
ylabel(YLabel,'fontsize',14);
hold on
plot_a = plot(t,VA,'-','LineWidth',2.5);
plot_b = plot(t,VB,'-','LineWidth',2.5);
set(plot_a,'Color',color_a);
set(plot_b,'Color',color_b);
axes = gca;
axes.LineWidth = 1.0;
axes.GridAlpha = 0.4;
grid on
hold off

% Compute an approximation for the first derivatives of the optimization
% variables of the model with respect to the varying exogenous parameter.
dVA = diff(VA);
dVB = diff(VB);
dt = diff(t);
dVAdt = dVA ./ dt;
dVBdt = dVB ./ dt;

% Set the name of the figure for the first derivatives.
FigureName = DerivativeFigureNameStrings{VariableIndex};
% Set the title of the figure for the first derivatives.
FigureTitle = sprintf('$$\\frac{d%s_A}{d%s}~\\textrm{vs}~\\frac{d%s_B}{d%s}$$',...
                      VariableName,ParamName,VariableName,ParamName);
% Set the x-axis label of the figure for the first derivatives.
XLabel = ParamName;
% Set the y-axis lable of the figure for the first derivatives.
YLabel = sprintf('$$\\frac{d%s_A}{d%s} / \\frac{d%s_B}{d%s}$$',...
                 VariableName,ParamName,VariableName,ParamName);
% Get the color for each curve for the first derivatives.
[color_a,color_b] = PlotColors(dVAdt,dVBdt,MinimumDigitsAccuracy);

% Do the actual plotting for the first derivative of the optimal
% optimization variables of the model with respect to the varying
% exogenous parameter.
figure('Name',FigureName);
title(FigureTitle,'fontsize',14,'Interpreter','latex');
xlabel(XLabel,'fontsize',14);
ylabel(YLabel,'fontsize',14,'Interpreter','latex');
hold on
plot_a = plot(t(2:end),dVAdt,'-','LineWidth',2.5);
plot_b = plot(t(2:end),dVBdt,'-','LineWidth',2.5);
set(plot_a,'Color',color_a);
set(plot_b,'Color',color_b);
axes = gca;
axes.LineWidth = 1.0;
axes.GridAlpha = 0.4;
grid on
hold off



end

function PlotEstimatedDerivatives(VariableIndex,ParamIndex,VA,VB,t,ExtraParams,MinimumDigitsAccuracy)

% Unfold the contents the exra parameters structure.
LA = ExtraParams.LA;
LB = ExtraParams.LB;
alpha = ExtraParams.alpha;
beta = ExtraParams.beta;
So = ExtraParams.So;
dTAdt = ExtraParams.dTAdt;
dTBdt = ExtraParams.dTBdt;

% Set the string names of the optimal endogenous parameters to be plotted 
% in Latex notation.
VariablesNames = {'T^{\ast}','s^{\ast}','x^{\ast}','p^{\ast}',...
                  'q^{\ast}','\Pi^{\ast}','R^{\ast}','C^{\ast}'};

% Set the string names of the exogenous parameters in Latex notation.
ParamsNames = {'\lambda_A','\lambda_B','p_A^{(0)}','p_B^{(0)}',...
               '\mu','\kappa','c','\gamma'};

% Set the estimated derivatives figure name strings.
EstimatedDerivativeFigureNameStrings = {'Optimal Investment Levels Estimated Derivatives',...
                                        'Optimal Limiting Influences Estimated Derivatives',...
                                        'Optimal Limiting Beliefs Estimated Derivatives',...
                                        'Optimal Prices Estimated Derivatives',...
                                        'Optimal Quantities Estimated Derivatives',...
                                        'Optimal Profits Estimated Derivatives',...
                                        'Optimal Revenues Estimated Derivatives',...
                                        'Optimal Costs Estimated Derivatives'};

% Set the name of the endogenous optimization variable to be plotted.
VariableName = VariablesNames{VariableIndex};
% Set the name of the exogenous model parameter to be plotted.
ParamName = ParamsNames{ParamIndex};
                                                                   
% Compute the estimated derivatives for the endogenous model variables:
% SA SB XA XB QA QB.

% Compute the estimated derivatives only for values of the VariableIndex
% in {2,3,5}.
switch VariableIndex
    % (i)  : theta_SA / theta_TA =  (LB/So) * (1 - SA)
    % (ii) : theta_SA / theta_TB = -(LA/So) * SA
    % (iii): theta_SB / theta_TB =  (LA/So) * (1 - SB)
    % (iv) : theta_SB / theat_TA = -(LB/So) * SB
    case 2
        theta_SA_theta_TA =  (LB ./ So) .* (1 - VA);
        theta_SA_theta_TB = -(LA ./ So) .* VA;
        theta_SB_theta_TB =  (LA ./ So) .* (1 - VB);
        theta_SB_theta_TA = -(LB ./ So) .* VB;
        dSAdt = theta_SA_theta_TA(2:end) .* dTAdt + theta_SA_theta_TB(2:end) .* dTBdt;
        dSBdt = theta_SB_theta_TA(2:end) .* dTAdt + theta_SB_theta_TB(2:end) .* dTBdt;
        
        % Set the name of the figure for the estimated first derivatives.
        FigureName = EstimatedDerivativeFigureNameStrings{VariableIndex};
        % Set the title of the figure for the estimated first derivatives.
        FigureTitle = sprintf('$$\\frac{d%s_A}{d%s}~\\textrm{vs}~\\frac{d%s_B}{d%s}$$',...
                      VariableName,ParamName,VariableName,ParamName);
        % Set the x-axis label of the figure for the estimated first derivatives.
        XLabel = ParamName;
        % Set the y-axis lable of the figure for the estimated first derivatives.
        YLabel = sprintf('$$\\frac{d%s_A}{d%s} / \\frac{d%s_B}{d%s}$$',...
                 VariableName,ParamName,VariableName,ParamName);
        % Get the color for each curve for the estimated first derivatives.
        [color_a,color_b] = PlotColors(dSAdt,dSBdt,MinimumDigitsAccuracy);

        % Do the actual plotting for the estimated first derivative of the optimal
        % optimization variables of the model with respect to the varying
        % exogenous parameter.
        figure('Name',FigureName);
        title(FigureTitle,'fontsize',14,'Interpreter','latex');
        xlabel(XLabel,'fontsize',14);
        ylabel(YLabel,'fontsize',14,'Interpreter','latex');
        hold on
        plot_a = plot(t(1:end-1),dSAdt,'-','LineWidth',2.5);
        plot_b = plot(t(1:end-1),dSBdt,'-','LineWidth',2.5);
        set(plot_a,'Color',color_a);
        set(plot_b,'Color',color_b);
        axes = gca;
        axes.LineWidth = 1.0;
        axes.GridAlpha = 0.4;
        grid on
        hold off
    % (i)  : theta_XA / theta_TA =  (LB/So) * (1 - XA)
    % (ii) : theta_XA / theta_TB = -(LA/So) * XA
    % (iii): theta_XB / theta_TB =  (LA/So) * (1 - XB)
    % (iv) : theta_XB / theat_TA = -(LB/So) * XB
    case 3
        theta_XA_theta_TA =  (LB ./ So) .* (1 - VA);
        theta_XA_theta_TB = -(LA ./ So) .* VA;
        theta_XB_theta_TB =  (LA ./ So) .* (1 - VB);
        theta_XB_theta_TA = -(LB ./ So) .* VB;
        dXAdt = theta_XA_theta_TA(2:end) .* dTAdt + theta_XA_theta_TB(2:end) .* dTBdt;
        dXBdt = theta_XB_theta_TA(2:end) .* dTAdt + theta_XB_theta_TB(2:end) .* dTBdt;
        
        % Set the name of the figure for the estimated first derivatives.
        FigureName = EstimatedDerivativeFigureNameStrings{VariableIndex};
        % Set the title of the figure for the estimated first derivatives.
        FigureTitle = sprintf('$$\\frac{d%s_A}{d%s}~\\textrm{vs}~\\frac{d%s_B}{d%s}$$',...
                      VariableName,ParamName,VariableName,ParamName);
        % Set the x-axis label of the figure for the estimated first derivatives.
        XLabel = ParamName;
        % Set the y-axis lable of the figure for the estimated first derivatives.
        YLabel = sprintf('$$\\frac{d%s_A}{d%s} / \\frac{d%s_B}{d%s}$$',...
                 VariableName,ParamName,VariableName,ParamName);
        % Get the color for each curve for the estimated first derivatives.
        [color_a,color_b] = PlotColors(dXAdt,dXBdt,MinimumDigitsAccuracy);

        % Do the actual plotting for the estimated first derivative of the optimal
        % optimization variables of the model with respect to the varying
        % exogenous parameter.
        figure('Name',FigureName);
        title(FigureTitle,'fontsize',14,'Interpreter','latex');
        xlabel(XLabel,'fontsize',14);
        ylabel(YLabel,'fontsize',14,'Interpreter','latex');
        hold on
        plot_a = plot(t(1:end-1),dXAdt,'-','LineWidth',2.5);
        plot_b = plot(t(1:end-1),dXBdt,'-','LineWidth',2.5);
        set(plot_a,'Color',color_a);
        set(plot_b,'Color',color_b);
        axes = gca;
        axes.LineWidth = 1.0;
        axes.GridAlpha = 0.4;
        grid on
        hold off
end
end