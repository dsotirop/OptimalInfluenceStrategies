function plot_parameters_tuples(Topts,Sopts,Xopts,Popts,Qopts,Fopts,Params,ParamIndex,ConstIndices,ConstValues)

% This function provides fundamental 2-D plotting operations concerning  
% the main output quantities of the oligopolistic optimal influences 
% model with respect to the selected varying parameter.

% The quantities of interest to be plotted are the following:
% [i]:   Optimal investement levels (Topts = [T1_A_opt,T2_A_opt,T1_B_opt,T2_B_opt]).
% [ii]:  Optimal limiting influences (Sopts = [SA_opt,S1_opt,S2_opt,SB_opt]).
% [iii]: Optimal limiting beliefs (Xopts  = [XA_opt,XB_opt]).
% [iv]:  Optimal prices (Popts = [pA_opt,pB_opt]).
% [v]:   Optimal quantities (Qopts = [Q_A_opt,Q_B_opt]).
% [vi]:  Optimal profits (Fopts = [F_A_opt,F_B_opt]).

% The fundamental plotting axes may correspond to the following external
% optimization parameters:
% [1]:  P_A_1 (ParamIndex = 1).
% [2]:  P_A_2 (ParamIndex = 2).
% [3]:  P_B_1 (ParamIndex = 3).
% [4]:  P_B_2 (ParamIndex = 4).
% [5]:  Lambda_A_1 (ParamIndex = 5).
% [6]:  Lambda_A_2 (ParamIndex = 6).
% [7]:  Lambda_B_1 (ParamIndex = 7).
% [8]:  Lambda_B_2 (ParamIndex = 8).
% [9]:  Theta1 (ParamIndex = 9).
% [10]: Theta2 (ParamIndex = 10).
% [11]: M (ParamIndex = 11).
% [12]: K (ParamIndex = 12).
% [13]: C (ParamIndex = 13).
% [14]: Gamma (ParamIndex = 14).

% [ParamIndex]: is the index of the parameter which is allowed to vary taking 
%               values within the discrete interval {1,...,14} .
% [Params]: is the matrix which stores the different 14-tuples in a row-wise 
%           manner for a given experimentation scenario. Mind that, only 
%           one of the columns will be allowed to vary.
% [ConstIndices]: is a vector storing the indices of the parameters which
%                 will be remaining constant. Therefore, it will hold that:
%                 ConstIndices = {1,..,14} \ {ParamIndex}.
% [ConstValues]: is a vector which stores the corresponding value for each 
%                one of the external optimization parameters.

% Construct a cell array storing the names of model parameters.
ParamsNames = {'P_{A,1}','P_{A,2}','P_{B,1}','P_{B,2}','Lambda_{A,1}','Lambda_{A,2}',...
               'Lambda_{B,1}','Lambda_{B,2}','Theta_1','Theta_2','M','K','C','Gamma'};
           
% Get the optimal investment levels for each different consumer-firm pair.          
T1_A_opt = Topts(:,1);
T2_A_opt = Topts(:,2);
T1_B_opt = Topts(:,3);
T2_B_opt = Topts(:,4);
% Get the limiting ifluence values for each agent {consumer or firm} in the network.
SA_opt = Sopts(:,1);
S1_opt = Sopts(:,2);
S2_opt = Sopts(:,3);
SB_opt = Sopts(:,4);
% Get the limiting beliefs for each product.
XA_opt = Xopts(:,1);
XB_opt = Xopts(:,2);
% Get the optimal prices for each product.
pA_opt = Popts(:,1);
pB_opt = Popts(:,2);
% Get the optimal quantities for each product.
Q_A_opt = Qopts(:,1);
Q_B_opt = Qopts(:,2);
% Get the optimal profits for each firm.
F_A_opt = Fopts(:,1);
F_B_opt = Fopts(:,2);
% Get the optimal revenues and costs for each firm.
F_A_Rev_opt = Fopts(:,3);
F_A_Cost_opt = Fopts(:,4);
F_B_Rev_opt = Fopts(:,5);
F_B_Cost_opt = Fopts(:,6);

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
    end;
    TitleString = strcat([TitleString ParamsNames{const_index} ' = ' num2str(ConstValues(const_index)) '|']);
end;
TitleName = sprintf(TitleString);

% Store the values and name of the parameter that is allowed to vary.
Param = Params(:,ParamIndex);
ParamName = ParamsNames{ParamIndex};

% 1st Figure: Plot T1_A_opt and T2_A_opt with respect to the varying parameter. 
Figure1Name  = strcat(['T1_A_opt and T2_A_opt Optimal with respect to ' ParamName]);
figure('Name',Figure1Name);
hold on
plot(Param,T1_A_opt,'*-r','LineWidth',1.4);
plot(Param,T2_A_opt,'*-g','LineWidth',1.4);
xlabel(ParamName);
ylabel('T1A_{opt} / T2A_{opt}');
grid on
title(TitleName);
% 2nd Figure: Plot T1_B_opt and T2_B_opt with respect to the varying parameter. 
Figure2Name  = strcat(['T1_B_opt and T2_B_opt Optimal with respect to ' ParamName]);
figure('Name',Figure2Name);
hold on
plot(Param,T1_B_opt,'*-r','LineWidth',1.4);
plot(Param,T2_B_opt,'*-g','LineWidth',1.4);
xlabel(ParamName);
ylabel('T1B_{opt} / T2B_{opt}');
grid on
title(TitleName);
% 3rd Figure: Plot SA_opt, S1_opt, S2_opt and SB_opt  with respect to the varying parameter.
Figure3Name  = strcat(['SA_opt, S1_opt, S2_opt and SB_opt Optimal with respect to ' ParamName]);
figure('Name',Figure3Name);
hold on
plot(Param,SA_opt,'*-c','LineWidth',1.4);
plot(Param,S1_opt,'*-r','LineWidth',1.4);
plot(Param,S2_opt,'*-g','LineWidth',1.4);
plot(Param,SB_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('SA_{opt} / S1_{opt} / S2_{opt} / SB_{opt}');
grid on
title(TitleName);
% 4th Figure: Plot XA_opt and XB_opt with respect to the varying parameter.
Figure4Name  = strcat(['XA_opt and XB_opt Optimal with respect to ' ParamName]);
figure('Name',Figure4Name);
subplot(2,1,1)
plot(Param,XA_opt,'*-c','LineWidth',1.4);
xlabel(ParamName);
ylabel('XA_{opt}');
grid on
title(TitleName);
subplot(2,1,2)
plot(Param,XB_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('XB_{opt}');
grid on
% 5th Figure: Plot pA_opt and pB_opt with respect to the varying parameter.
Figure5Name  = strcat(['pA_opt and pB_opt Optimal with respect to ' ParamName]);
figure('Name',Figure5Name);
subplot(2,1,1)
plot(Param,pA_opt,'*-c','LineWidth',1.4);
xlabel(ParamName);
ylabel('pA_{opt}');
grid on
title(TitleName);
subplot(2,1,2)
plot(Param,pB_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('pB_{opt}');
grid on
% 6th Figure: Plot Q_A_opt and Q_B_opt with respect to the varying parameter.
Figure6Name  = strcat(['Q_A_opt and Q_B_opt Optimal with respect to ' ParamName]);
figure('Name',Figure6Name);
subplot(2,1,1)
plot(Param,Q_A_opt,'*-c','LineWidth',1.4);
xlabel(ParamName);
ylabel('QA_{opt}');
grid on
title(TitleName);
subplot(2,1,2)
plot(Param,Q_B_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('QB_{opt}');
grid on
% 7th Figure: Plot F_A_opt and F_B_opt with respect to the varying parameter.
Figure7Name  = strcat(['F_A_opt and F_B_opt Optimal with respect to ' ParamName]);
figure('Name',Figure7Name);
subplot(2,1,1)
plot(Param,F_A_opt,'*-c','LineWidth',1.4);
xlabel(ParamName);
ylabel('FA_{opt}');
grid on
title(TitleName);
subplot(2,1,2)
plot(Param,F_B_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('FB_{opt}');
grid on
% 8th Figure: Plot F_A_Rev_opt and F_A_Cost_opt with respect to the varying parameter.
Figure8Name  = strcat(['F_A_Rev_opt and F_A_Cost_opt Optimal with respect to ' ParamName]);
figure('Name',Figure8Name);
hold on
plot(Param,F_A_Rev_opt,'*-b','LineWidth',1.4);
plot(Param,F_A_Cost_opt,'*-r','LineWidth',1.4);
plot(Param,F_A_opt,'*-c','LineWidth',1.4);
xlabel(ParamName);
ylabel('FA-Rev_{opt} / FA-Cost_{opt}');
grid on
title(TitleName);
% 9th Figure: Plot F_B_Rev_opt and F_B_Cost_opt with respect to the varying parameter.
Figure8Name  = strcat(['F_B_Rev_opt and F_B_Cost_opt Optimal with respect to ' ParamName]);
figure('Name',Figure8Name);
hold on
plot(Param,F_B_Rev_opt,'*-b','LineWidth',1.4);
plot(Param,F_B_Cost_opt,'*-r','LineWidth',1.4);
plot(Param,F_B_opt,'*-m','LineWidth',1.4);
xlabel(ParamName);
ylabel('FB-Rev_{opt} / FB-Cost_{opt}');
grid on
title(TitleName);
end

