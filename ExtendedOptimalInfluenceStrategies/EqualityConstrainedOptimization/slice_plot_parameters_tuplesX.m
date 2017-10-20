function slice_plot_parameters_tuplesX(TLopts,Fvals,Sopts,X,Params,ParamIndices,ParamSubRanges,ConstIndices,ConstValues)

% This function performs a series of 2-dimensional plotting operations by
% slicing the two-dimensional grid of paramaters stored in the submatrix
% Params(ParamIndices,:) in such a way so that one of the parameters will
% remain constant at each plot taking values that are determined by the contents
% of the matrix ParamSubRanges. The input argument ParamSubRanges is
% assumed to be a 2-rowed matrix of the following form:
%                  | Param1_1 ... Param1_N |
% ParamSubRanges = |
%                  | Param2_1 ... Param2_N |
% This entails that the maximum number of discrete parameter values for the
% the pair of parameters under investigation will be N.
% Mind that the number of varying parameters must be necessarily 2 for this
% function to work properly.

% Get the number of discrete parameter values for each parameter.
N = size(ParamSubRanges,2);

% Store input arguements into separate variables.
Param1Values = ParamSubRanges(1,:);
Param2Values = ParamSubRanges(2,:);

% Construct a cell array storing the names of model parameters.
ParamsNames = {'P0(0)','P1(0)','P2(0)','Theta1','Theta2','Delta','Gamma','K','Beta'};

% Get the optimal values for T1,T2,Lambda1 and Lambda2.
T1 = TLopts(:,1);
T2 = TLopts(:,2);
Lambda1 = TLopts(:,3);
Lambda2 = TLopts(:,4);
% Get the optimal values for S0, S1 and S2.
S0 = Sopts(:,1);
S1 = Sopts(:,2);
S2 = Sopts(:,3);

% Store separately the indices of the 7 parameters that remain constant.
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

% Set the marker specifiers.
MarkerSpecifiers = {'*','.','x','+','o','s','d'};

% If the number of discrete paramater values is less than or equal to 7
% then all the 2d plotting operations for each free parameter with respect
% to other (which will be remaining constant) will be appearing on the same
% figure.
if(N<=7)
    % 1st Figure.
    Figure1Name = strcat(['T1 / T2 Optimal with respect to ' Param1Name]);
    figure('Name',Figure1Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param2Indices = find(Param2==Param2Values(n));
        param1_values = Param1(Param2Indices);
        T1_values = T1(Param2Indices);
        T2_values = T2(Param2Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        plot(param1_values,T1_values,line1_spec,'LineWidth',1.4);
        plot(param1_values,T2_values,line2_spec,'LineWidth',1.4);
        legend_string = strcat([Param2Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string];
    end;
    xlabel(Param1Name);
    ylabel('T1opt / T2opt');
    legend(LegendString);
    grid on
    hold off
    title(TitleName);
    % 2nd Figure.
    Figure2Name = strcat(['T1 / T2 Optimal with respect to ' Param2Name]);
    figure('Name',Figure2Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param1Indices = find(Param1==Param1Values(n));
        param2_values = Param2(Param1Indices);
        T1_values = T1(Param1Indices);
        T2_values = T2(Param1Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        plot(param2_values,T1_values,line1_spec,'LineWidth',1.4);
        plot(param2_values,T2_values,line2_spec,'LineWidth',1.4);
        legend_string = strcat([Param1Name ' = ' num2str(Param1Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string];
    end;
    xlabel(Param2Name);
    ylabel('T1opt / T2opt');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    % 3rd Figure.
    Figure3Name = strcat(['S0 / S1 / S2 Optimal with respect to ' Param1Name]);
    figure('Name',Figure3Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param2Indices = find(Param2==Param2Values(n));
        param1_values = Param1(Param2Indices);
        S0_values = S0(Param2Indices);
        S1_values = S1(Param2Indices);
        S2_values = S2(Param2Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        line3_spec = strcat([mark_spec 'b']);
        plot(param1_values,S1_values,line1_spec,'LineWidth',1.4);
        plot(param1_values,S2_values,line2_spec,'LineWidth',1.4);
        plot(param1_values,S0_values,line3_spec,'LineWidth',1.4);
        legend_string = strcat([Param2Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string;legend_string];
    end;
    xlabel(Param1Name);
    ylabel('S0opt / S1opt / S2opt');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    %4th Figure.
    Figure4Name = strcat(['S0 / S1 / S2 Optimal with respect to ' Param2Name]);
    figure('Name',Figure4Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param1Indices = find(Param1==Param1Values(n));
        param2_values = Param2(Param1Indices);
        S0_values = S0(Param1Indices);
        S1_values = S1(Param1Indices);
        S2_values = S2(Param1Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        line3_spec = strcat([mark_spec 'b']);
        plot(param2_values,S1_values,line1_spec,'LineWidth',1.4);
        plot(param2_values,S2_values,line2_spec,'LineWidth',1.4);
        plot(param2_values,S0_values,line3_spec,'LineWidth',1.4);
        legend_string = strcat([Param1Name ' = ' num2str(Param1Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string;legend_string];
    end;
    xlabel(Param2Name);
    ylabel('S0opt / S1opt / S2opt');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    % 5th Figure.
    Figure5Name = strcat(['Limiting Beliefs with respect to ' Param1Name]);
    figure('Name',Figure5Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param2Indices = find(Param2==Param2Values(n));
        param1_values = Param1(Param2Indices);
        X_values = X(Param2Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'm']);
        plot(param1_values,X_values,line1_spec,'LineWidth',1.4);
        legend_string = strcat([Param2Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string];
    end;
    xlabel(Param1Name);
    ylabel('Limiting Beliefs');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    % 6th Figure.
    Figure6Name = strcat(['Limiting Beliefs with respect to ' Param2Name]);
    figure('Name',Figure6Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param1Indices = find(Param1==Param1Values(n));
        param2_values = Param2(Param1Indices);
        X_values = X(Param1Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'm']);
        plot(param2_values,X_values,line1_spec,'LineWidth',1.4);
        legend_string = strcat([Param1Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string];
    end;
    xlabel(Param2Name);
    ylabel('Limiting Beliefs');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    % 7th Figure
    Figure7Name = strcat(['Optimal Profit with respect to ' Param1Name]);
    figure('Name',Figure7Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param2Indices = find(Param2==Param2Values(n));
        param1_values = Param1(Param2Indices);
        F_values = Fvals(Param2Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'b']);
        plot(param1_values,F_values,line1_spec,'LineWidth',1.4);
        legend_string = strcat([Param2Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string];
    end;
    xlabel(Param1Name);
    ylabel('Optimal Profit');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off
    % 8th Figure
    Figure8Name = strcat(['Optimal Profit with respect to ' Param2Name]);
    figure('Name',Figure8Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param1Indices = find(Param1==Param1Values(n));
        param2_values = Param2(Param1Indices);
        F_values = Fvals(Param1Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'b']);
        plot(param2_values,F_values,line1_spec,'LineWidth',1.4);
        legend_string = strcat([Param1Name ' = ' num2str(Param1Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string];
    end;
    xlabel(Param2Name);
    ylabel('Optimal Profit');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off

    % 9th Figure.
    Figure9Name = strcat(['Lambda1 / Lambda2 Optimal with respect to ' Param1Name]);
    figure('Name',Figure9Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param2Indices = find(Param2==Param2Values(n));
        param1_values = Param1(Param2Indices);
        Lambda1_values = Lambda1(Param2Indices);
        Lambda2_values = Lambda2(Param2Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        plot(param1_values,Lambda1_values,line1_spec,'LineWidth',1.4);
        plot(param1_values,Lambda2_values,line2_spec,'LineWidth',1.4);
        legend_string = strcat([Param2Name ' = ' num2str(Param2Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string];
    end;
    xlabel(Param1Name);
    ylabel('Lambda1opt / Lambda2opt');
    legend(LegendString);
    grid on
    hold off
    title(TitleName);
    % 10th Figure.
    Figure10Name = strcat(['Lambda1 / Lambda2 Optimal with respect to ' Param2Name]);
    figure('Name',Figure10Name);
    LegendString = [];
    hold on
    for n = 1:1:N
        Param1Indices = find(Param1==Param1Values(n));
        param2_values = Param2(Param1Indices);
        Lambda1_values = Lambda1(Param1Indices);
        Lambda2_values = Lambda2(Param1Indices);
        mark_spec = MarkerSpecifiers{n};
        line1_spec = strcat([mark_spec 'r']);
        line2_spec = strcat([mark_spec 'g']);
        plot(param2_values,Lambda1_values,line1_spec,'LineWidth',1.4);
        plot(param2_values,Lambda2_values,line2_spec,'LineWidth',1.4);
        legend_string = strcat([Param1Name ' = ' num2str(Param1Values(n),'%.2f')]);
        LegendString = [LegendString;legend_string;legend_string];
    end;
    xlabel(Param2Name);
    ylabel('Lambda1opt / Lambda2opt');
    legend(LegendString);
    grid on
    title(TitleName);
    hold off    
end;        


end

