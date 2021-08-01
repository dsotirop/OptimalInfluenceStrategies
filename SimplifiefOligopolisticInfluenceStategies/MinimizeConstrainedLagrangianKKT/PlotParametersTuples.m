function PlotParametersTuples(Topts,Sopts,Xopts,Popts,Qopts,Fopts,Params,ParamIndex,ConstIndices,ConstValues,MinimumDigitsAccuracy)

% This function provides fundamental 2-D plotting operations concerning  
% the main output quantities of the simplified oligopolistic optimal  
% influences model with respect to the selected varying parameter.

% The quantities of interest to be plotted are the following:
% [i]:   Optimal investement levels (Topts = [TAopt,TBopt]).
% [ii]:  Optimal limiting influences (Sopts = [SAopt,SCopt,SBopt]).
% [iii]: Optimal limiting beliefs (Xopts  = [XAopt,XBopt]).
% [iv]:  Optimal prices (Popts = [PAopt,PBopt]).
% [v]:   Optimal quantities (Qopts = [QAopt,QBopt]).
% [vi]:  Optimal profits (Fopts = [FAopt,FBopt,FA_rev_opt,FB_rev_opt,FA_cost_opt,FB_cost_opt]).

% The fundamental plotting axes may correspond to the following external
% optimization parameters:
% [1]: LA (ParamIndex = 1).
% [2]: LB (ParamIndex = 2).
% [3]: PA (ParamIndex = 3).
% [4]: PB (ParamIndex = 4).
% [5]: M  (ParamIndex = 5). 
% [6]: K  (ParamIndex = 6).
% [7]: C  (ParamIndex = 7).
% [8]: G  (ParamIndex = 8).

% [ParamIndex]: is the index of the parameter which is allowed to vary taking 
%               values within the discrete interval {1,...,8} .
% [Params]: is the matrix which stores the different 8-tuples in a row-wise 
%           manner for a given experimentation scenario. Mind that, only 
%           one of the columns will be allowed to vary.
% [ConstIndices]: is a vector storing the indices of the parameters which
%                 will be remaining constant. Therefore, it will hold that:
%                 ConstIndices = {1,..,8} \ {ParamIndex}.
% [ConstValues]: is a vector which stores the corresponding value for each 
%                one of the external optimization parameters.

% Construct a cell array storing the names of model parameters.
ParamsNames = {'LA','LB','PA','PB','M','K','C','G'};

% Get the optimal investment levels for each firm.
TAopt = Topts(:,1);
TBopt = Topts(:,2);
% Get the valid indices for which the fundamental optimization variables is
% different than -1. The current code implementation assumes that when the
% underlying optimization problem admits no valid solutions then the main 
% optimization variables will be assigned the values: TAopt=-1 and TBop=-1.
Ivalid = (TAopt~=-1);
TAopt = TAopt(Ivalid);
TBopt = TBopt(Ivalid);
% Get the limiting ifluence values for each agent {consumer or firm} in the network.
SAopt = Sopts(:,1);
SCopt = Sopts(:,2);
SBopt = Sopts(:,3);
SAopt = SAopt(Ivalid);
SCopt = SCopt(Ivalid);
SBopt = SBopt(Ivalid);
% Get the limiting beliefs for each product.
XAopt = Xopts(:,1);
XBopt = Xopts(:,2);
XAopt = XAopt(Ivalid);
XBopt = XBopt(Ivalid);
% Get the optimal prices for each product.
PAopt = Popts(:,1);
PBopt = Popts(:,2);
PAopt = PAopt(Ivalid);
PBopt = PBopt(Ivalid);
% Get the optimal quantities for each product.
QAopt = Qopts(:,1);
QBopt = Qopts(:,2);
QAopt = QAopt(Ivalid);
QBopt = QBopt(Ivalid);
% Get the optimal profits for each firm.
FAopt = Fopts(:,1);
FBopt = Fopts(:,2);
FAopt = FAopt(Ivalid);
FBopt = FBopt(Ivalid);
% Get the optimal revenues and costs for each firm.
FA_rev_opt = Fopts(:,3);
FB_rev_opt = Fopts(:,4);
FA_cost_opt = Fopts(:,5);
FB_cost_opt = Fopts(:,6);
FA_rev_opt = FA_rev_opt(Ivalid);
FB_rev_opt = FB_rev_opt(Ivalid);
FA_cost_opt = FA_cost_opt(Ivalid);
FB_cost_opt = FB_cost_opt(Ivalid);
% Compute the competition-related parameters alpha and beta.
M = Params(:,5);
K = Params(:,6);
ALPHA = (K.*M - 2) ./ (M.^2 - 4);
BETA = (2*K - M) ./ (M.^2 - 4);
ALPHA = ALPHA(Ivalid);
BETA = BETA(Ivalid);
% Set a variable for storing the number of constant parameters.
const_parameters_num = length(ConstIndices);
% Set the next line indicator. This variable will be storing the number of 
% const_parameter_name = const_parameter_value pairs that will be appearing 
% within the same line of the title string.
new_line_indicator = 4;

% Generate the title string for each figure.
% Initialize the title string.
TitleString = '';
for const_index = 1:1:const_parameters_num
    if(mod(const_index,new_line_indicator+1)==0)
        TitleString = strcat([TitleString '\n']);
    end
    TitleString = strcat([TitleString ParamsNames{ConstIndices(const_index)} ' = ' num2str(ConstValues(const_index)) '|']);
end
TitleName = sprintf(TitleString);

% Store the values and name of the parameter that is allowed to vary.
Params = Params(Ivalid,:);
Param = Params(:,ParamIndex);
ParamName = ParamsNames{ParamIndex};

% 1st Figure: Plot TAopt and TBopt with respect to the varying parameter. 
Figure1Name  = strcat(['TAopt and TBopt Optimal with respect to ' ParamName]);
figure('Name',Figure1Name);
hold on
plot(Param,TAopt,'.-r','LineWidth',1.8);
plot(Param,TBopt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('TA_{opt} / TB_{opt}');
grid on
title(TitleName);
% 2nd Figure: Plot SAopt, SCopt and SBopt with respect to the varying parameter.
Figure2Name  = strcat(['SAopt, SCopt, and SBopt Optimal with respect to ' ParamName]);
figure('Name',Figure2Name);
hold on
plot(Param,SAopt,'.-r','LineWidth',1.8);
plot(Param,SCopt,'.-b','LineWidth',1.8);
plot(Param,SBopt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('SA_{opt} / SC_{opt} / SB_{opt}');
grid on
title(TitleName);
% 3rd Figure: Plot XAopt and XBopt with respect to the varying parameter.
Figure3Name = strcat(['XAopt and XBopt Optimal with respect to ' ParamName]);
figure('Name',Figure3Name);
hold on
plot(Param,XAopt,'.-r','LineWidth',1.8);
plot(Param,XBopt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('XA_{opt} / XB_{opt}');
if(abs(min(XAopt)-max(XAopt))<10^(-MinimumDigitsAccuracy))
    axis([0 1 0 1]);
end
grid on
title(TitleName);
% 4th Figure: Plot PAopt and PBopt with respect to the varying parameter.
Figure4Name = strcat(['PAopt and PBopt Optimal with respect to ' ParamName]);
figure('Name',Figure4Name);
hold on
plot(Param,PAopt,'.-r','LineWidth',1.8);
plot(Param,PBopt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('PA_{opt} / PB_{opt}');
grid on
title(TitleName);
% 5th Figure: Plot QAopt and QBopt with respect to the varying parameter.
Figure5Name = strcat(['QAopt and QBopt Optimal with respect to ' ParamName]);
figure('Name',Figure5Name);
hold on
plot(Param,QAopt,'.-r','LineWidth',1.8);
plot(Param,QBopt,'.-g','LineWidth',1.8);
% plot(Param,ALPHA,'--m','LineWidth',1.8);
% plot(Param,BETA,'--b','LineWidth',1.8);
% plot(Param,ALPHA+BETA,'--k','LineWidth',1.8);
xlabel(ParamName);
ylabel('QA_{opt} / QB_{opt}');
grid on
title(TitleName);
% 6th Figure: Plot FAopt and FBopt with respect to the varying parameter.
Figure6Name = strcat(['FAopt and FBopt Optimal with respect to ' ParamName]);
figure('Name',Figure6Name);
hold on
plot(Param,FAopt,'.-r','LineWidth',1.8);
plot(Param,FBopt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('FA_{opt} / FB_{opt}');
grid on
title(TitleName);
% 7th Figure: Plot FA_rev_opt and FB_rev_opt with respect to the varying parameter.
Figure7Name = strcat(['FA_rev_opt and FB_rev_opt Optimal with respect to ' ParamName]);
figure('Name',Figure7Name);
hold on
plot(Param,FA_rev_opt,'.-r','LineWidth',1.8);
plot(Param,FB_rev_opt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('FA_{rev}_{opt} / FB_{rev}_{opt}');
grid on
title(TitleName);
% 8th Figure: Plot FA_cost_opt and FB_cost_opt with respect to the varying parameter.
Figure8Name = strcat(['FA_cost_opt and FB_cost_opt Optimal with respect to ' ParamName]);
figure('Name',Figure8Name);
hold on
plot(Param,FA_cost_opt,'.-r','LineWidth',1.8);
plot(Param,FB_cost_opt,'.-g','LineWidth',1.8);
xlabel(ParamName);
ylabel('FA_{cost}_{opt} / FB_{cost}_{opt}');
grid on
title(TitleName);
end