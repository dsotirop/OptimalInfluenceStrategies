function plot_parameter_triples(Topts,Fvals,Params,ParamIndices,ConstIndices,ConstValues)

% This function plots in 3 dimensions the optimal T1 and T2 with respect
% to the parameters indexed by ParamIndex1 and ParamIndex2.

ParamIndex1 = ParamIndices(1);
ParamIndex2 = ParamIndices(2);
ConstIndex1 = ConstIndices(1);
ConstIndex2 = ConstIndices(2);
ConstIndex3 = ConstIndices(3);

ParamsNames = {'P1(0)','P2(0)','Theta','Delta','Gamma'};
T1 = Topts(:,1);
T2 = Topts(:,2);

Param1 = Params(:,ParamIndex1);
Param1Name = ParamsNames{ParamIndex1};
Param2 = Params(:,ParamIndex2);
Param2Name = ParamsNames{ParamIndex2};

Const1Name = ParamsNames{ConstIndex1};
Const2Name = ParamsNames{ConstIndex2};
Const3Name = ParamsNames{ConstIndex3};

Const1ValueString = num2str(ConstValues(1));
Const2ValueString = num2str(ConstValues(2));
Const3ValueString = num2str(ConstValues(3));

Figure1Name = strcat(['T1 Optimal with respect to ' Param1Name ' and ' Param2Name]);
TitleName = strcat([Const1Name ' = ' Const1ValueString ' ' Const2Name ' = ' Const2ValueString ' ' Const3Name ' = ' Const3ValueString ' ']);
figure('Name',Figure1Name);
plot3(Param1,Param2,T1,'*r','LineWidth',1.4);
xlabel(Param1Name);
ylabel(Param2Name);
zlabel('T1opt');
grid on
title(TitleName);

Figure2Name = strcat(['T2 Optimal with respect to ' Param1Name ' and ' Param2Name]);
figure('Name',Figure2Name);
plot3(Param1,Param2,T2,'*g','LineWidth',1.4);
xlabel(Param1Name);
ylabel(Param2Name);
zlabel('T2opt');
grid on
title(TitleName);

Figure3Name = strcat(['Optimal Profit with respect to ' Param1Name ' and ' Param2Name]);
figure('Name',Figure3Name);
plot3(Param1,Param2,Fvals,'*b','LineWidth',1.4);
xlabel(Param1Name);
ylabel(Param2Name);
zlabel('Profit');
grid on
title(TitleName);

Figure4Name = strcat(['Optimal T1 and T2 with respect to ' Param1Name ' and ' Param2Name]);
figure('Name',Figure3Name);
plot3(Param1,Param2,T1,'*r','LineWidth',1.4);
hold on
plot3(Param1,Param2,T2,'*g','LineWidth',1.4);
xlabel(Param1Name);
ylabel(Param2Name);
zlabel('T1opt / T2opt');
grid on
title(TitleName);

end