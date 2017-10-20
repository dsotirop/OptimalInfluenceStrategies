function plot_parameters_tuples(Topts,Fvals,Frevenue,Fcost,Sopts,X,Params,ParamIndices,ConstIndices,ConstValues)

% This function performs 2-dimensional or 3-dimensional plotting operations
% based on the number of free parameters that are indicated by the contents
% of the ParamIndices vector. The quantities of interest correspond to the
% input arguments Topts (optimal connection strengths), Fvals (optimal profit),
% Sopts (optimal limiting influences) and X (consensus belief).

% Construct a cell array storing the names of model parameters.
ParamsNames = {'P0(0)','P1(0)','P2(0)','Lambda1','Lambda2','Theta1','Theta2','Delta','Gamma'};
% Get the number of model parameters that are allowed to vary.
varying_parameters_num = length(ParamIndices);
% Get the optimal values for T1 and T2.
T1 = Topts(:,1);
T2 = Topts(:,2);
% Get the optimal values for S0, S1 and S2.
S0 = Sopts(:,1);
S1 = Sopts(:,2);
S2 = Sopts(:,3);
% Perform 2-dimensional or 3-dimensional plotting depending on the number
% of free parameters.
if(varying_parameters_num==1)
    % Store separately the indices of the 8 parameters that remain constant.
    ConstIndex1 = ConstIndices(1);
    ConstIndex2 = ConstIndices(2);
    ConstIndex3 = ConstIndices(3);
    ConstIndex4 = ConstIndices(4);
    ConstIndex5 = ConstIndices(5);
    ConstIndex6 = ConstIndices(6);
    ConstIndex7 = ConstIndices(7);
    ConstIndex8 = ConstIndices(8);
    % Store separately the names of the 8 parameters that remain constant.
    Const1Name = ParamsNames{ConstIndex1};
    Const2Name = ParamsNames{ConstIndex2};
    Const3Name = ParamsNames{ConstIndex3};
    Const4Name = ParamsNames{ConstIndex4};
    Const5Name = ParamsNames{ConstIndex5};
    Const6Name = ParamsNames{ConstIndex6};
    Const7Name = ParamsNames{ConstIndex7};
    Const8Name = ParamsNames{ConstIndex8};
    % Store separately the value strings of the 8 parameters that remain constant.
    Const1ValueString = num2str(ConstValues(1));
    Const2ValueString = num2str(ConstValues(2));
    Const3ValueString = num2str(ConstValues(3));
    Const4ValueString = num2str(ConstValues(4));
    Const5ValueString = num2str(ConstValues(5));
    Const6ValueString = num2str(ConstValues(6));
    Const7ValueString = num2str(ConstValues(7));
    Const8ValueString = num2str(ConstValues(8));
    % Store the values and name of the parameter that is allowed to vary.
    ParamIndex = ParamIndices(1);
    Param = Params(:,ParamIndex);
    ParamName = ParamsNames{ParamIndex};
    % Generate the title string for each figure.
    TitleString = strcat([Const1Name ' = ' Const1ValueString ' ' Const2Name ' = ' ...
                        Const2ValueString ' ' Const3Name ' = ' Const3ValueString ' ' ...
                        Const4Name ' = ' Const4ValueString '\n' Const5Name ' = ' Const5ValueString ' '...
                        Const6Name ' = ' Const6ValueString ' ' Const7Name ' = ' Const7ValueString ' '...
                        Const8Name ' = ' Const8ValueString]);
    TitleName = sprintf(TitleString);
    % 1st Figure.
    Figure1Name = strcat(['T1 and T2 Optimal with respect to ' ParamName]);
    figure('Name',Figure1Name);
    hold on
    plot(Param,T1,'*r','LineWidth',1.4);
    plot(Param,T2,'*g','LineWidth',1.4);
    xlabel(ParamName);
    ylabel('T1opt / T2opt');
    grid on
    title(TitleName);
    % 2nd Figure.
    Figure2Name = strcat(['Optimal Profit with respect to ' ParamName]);
    figure('Name',Figure2Name);
    plot(Param,Fvals,'*b','LineWidth',1.4);
    xlabel(ParamName);
    ylabel('Profit');
    grid on
    title(TitleName);
    % 3rd Figure.
    Figure3Name = strcat(['Optimal Influences with respect to ' ParamName]);
    figure('Name',Figure3Name);
    hold on
    plot(Param,S0,'*b','LineWidth',1.4);
    plot(Param,S1,'*r','LineWidth',1.4);
    plot(Param,S2,'*g','LineWidth',1.4);
    xlabel(ParamName);
    ylabel('S0opt / S1opt / S2opt');
    grid on
    title(TitleName);
    % 4th Figure.
    Figure4Name = strcat(['Limiting Beliefs with respect to ' ParamName]);
    figure('Name',Figure4Name);
    plot(Param,X,'*m','LineWidth',1.4);
    xlabel(ParamName);
    ylabel('Limiting Beliefs');
    grid on
    title(TitleName);
    % 5th Figure.
    Figure5Name = strcat(['Cost-Revenue Analysis with respect to ' ParamName]);
    figure('Name',Figure5Name);
    hold on
    plot(Param,Frevenue,'*m','LineWidth',1.4);
    plot(Param,Fvals,'*b','LineWidth',1.4);
    plot(Param,Fcost,'*k','LineWidth',1.4);
    xlabel(ParamName);
    ylabel('Revenue / Profit / Cost');
    grid on
    title(TitleName);
end;
if(varying_parameters_num==2)
    % Store separately the indices of the 6 parameters that remain constant.
    ConstIndex1 = ConstIndices(1);
    ConstIndex2 = ConstIndices(2);
    ConstIndex3 = ConstIndices(3);
    ConstIndex4 = ConstIndices(4);
    ConstIndex5 = ConstIndices(5);
    ConstIndex6 = ConstIndices(6);
    ConstIndex7 = ConstIndices(7);
    % Store separately the names of the 7 parameters that remain constant.
    Const1Name = ParamsNames{ConstIndex1};
    Const2Name = ParamsNames{ConstIndex2};
    Const3Name = ParamsNames{ConstIndex3};
    Const4Name = ParamsNames{ConstIndex4};
    Const5Name = ParamsNames{ConstIndex5};
    Const6Name = ParamsNames{ConstIndex6};
    Const7Name = ParamsNames{ConstIndex7};
    % Store separately the value strings of the 7 parameters that remain constant.
    Const1ValueString = num2str(ConstValues(1));
    Const2ValueString = num2str(ConstValues(2));
    Const3ValueString = num2str(ConstValues(3));
    Const4ValueString = num2str(ConstValues(4));
    Const5ValueString = num2str(ConstValues(5));
    Const6ValueString = num2str(ConstValues(6));
    Const7ValueString = num2str(ConstValues(7));
    % Store separately the values and names of the parameters that are
    % allowed to vary.
    ParamIndex1 = ParamIndices(1);
    ParamIndex2 = ParamIndices(2);
    Param1 = Params(:,ParamIndex1);
    Param1Name = ParamsNames{ParamIndex1};
    Param2 = Params(:,ParamIndex2);
    Param2Name = ParamsNames{ParamIndex2};
    % Generate the title string for each figure.
    TitleString = strcat([Const1Name ' = ' Const1ValueString ' ' Const2Name ' = ' ...
                        Const2ValueString ' ' Const3Name ' = ' Const3ValueString ' ' ...
                        Const4Name ' = ' Const4ValueString '\n' Const5Name ' = ' Const5ValueString ' '...
                        Const6Name ' = ' Const6ValueString ' ' Const7Name ' = ' Const7ValueString]);
    TitleName = sprintf(TitleString);
    % 1st Figure.
    Figure1Name = strcat(['T1 Optimal with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure1Name);
    plot3(Param1,Param2,T1,'*r','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('T1opt');
    grid on
    title(TitleName);
    % 2nd Figure.
    Figure2Name = strcat(['T2 Optimal with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure2Name);
    plot3(Param1,Param2,T2,'*g','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('T2opt');
    grid on
    title(TitleName);
    % 3rd Figure.
    Figure3Name = strcat(['T1 / T2 Optimal with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure3Name);
    plot3(Param1,Param2,T1,'*r','LineWidth',1.4);
    hold on
    plot3(Param1,Param2,T2,'*g','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('T1opt / T2opt');
    grid on
    title(TitleName);
    % 4th Figure.
    Figure4Name = strcat(['Optimal Profit with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure4Name);
    plot3(Param1,Param2,Fvals,'*b','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('Optimal Profit');
    grid on
    title(TitleName);
    % 5th Figure.
    Figure5Name = strcat(['S0 / S1 / S2 Optimal with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure5Name);
    plot3(Param1,Param2,S1,'*r','LineWidth',1.4);
    hold on
    plot3(Param1,Param2,S2,'*g','LineWidth',1.4);
    plot3(Param1,Param2,S0,'*b','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('S0opt / S1opt / S2opt');
    grid on
    title(TitleName);
    % 6th Figure.
    Figure6Name = strcat(['Limiting Beliefs with respect to ' Param1Name ' and ' Param2Name]);
    figure('Name',Figure6Name);
    plot3(Param1,Param2,X,'*m','LineWidth',1.4);
    xlabel(Param1Name);
    ylabel(Param2Name);
    zlabel('Limiting Beliefs');
    grid on
    title(TitleName);
end;
    
end

