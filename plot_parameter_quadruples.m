function plot_parameter_quadruples(Topts,Fvals,Params,ParamIndex,ConstIndices,ConstValues)

% This function plots in 3 dimensions the optimal T1 and T2 with respect
% to the parameters indexed by ParamIndex1 and ParamIndex2.

ConstIndex1 = ConstIndices(1);
ConstIndex2 = ConstIndices(2);
ConstIndex3 = ConstIndices(3);
ConstIndex4 = ConstIndices(4);

ParamsNames = {'P1(0)','P2(0)','Theta','Delta','Gamma'};
T1 = Topts(:,1);
T2 = Topts(:,2);

Param = Params(:,ParamIndex);
ParamName = ParamsNames{ParamIndex};

Const1Name = ParamsNames{ConstIndex1};
Const2Name = ParamsNames{ConstIndex2};
Const3Name = ParamsNames{ConstIndex3};
Const4Name = ParamsNames{ConstIndex4};

Const1ValueString = num2str(ConstValues(1));
Const2ValueString = num2str(ConstValues(2));
Const3ValueString = num2str(ConstValues(3));
Const4ValueString = num2str(ConstValues(4));

Figure1Name = strcat(['T1 and T2 Optimal with respect to ' ParamName]);
TitleName = strcat([Const1Name ' = ' Const1ValueString ' ' Const2Name ' = ' Const2ValueString ' ' Const3Name ' = ' Const3ValueString ' ' Const4Name ' = ' Const4ValueString]);
figure('Name',Figure1Name);
hold on
plot(Param,T1,'*r','LineWidth',1.4);
plot(Param,T2,'*g','LineWidth',1.4);
xlabel(ParamName);
ylabel('T1opt / T2opt');
grid on
title(TitleName);

Figure2Name = strcat(['Optimal Profit with respect to ' ParamName]);
figure('Name',Figure2Name);
plot(Param,Fvals,'*b','LineWidth',1.4);
xlabel(ParamName);
ylabel('Profit');
grid on
title(TitleName);


end

