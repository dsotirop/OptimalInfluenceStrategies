function plot_influence_regions(Topts,Sopts,Flags1,Flags2,Params,ParamIndices,ConstIndices,ConstValues,CommonParameterName)

% This function provides the fundamental functionality for plotting the
% influence regions (Theta_hat - Theta_star) as a function of one primary
% varying parameter other than Theta1 and Theta2 since influence regions
% can only be defined with respect to Theta parameters. The number of
% varying parameters defined in vector ParamIndices can be either 2 or 3.
% When the number of varying parameters is 3 there are two possible cases. 
% The 1st case corresponds to the situation where 2 of the varying parameters
% take exactly the same values. In this case they will be identified by the
% CommonParameterName. The 2nd case corresponds to the situation where one
% of the varying paramers takes values in a narrower interval than the
% second one. In this case a series of graphs will be produced for each
% value of the parameter that has the shorter range.
% -------------------------------------------------------------------------
% It is extremelly important to note that the index of the varying Theta
% parameter must always be stored within the last position of ParamIndices
% vector otherwise the function will not operate properly.
% -------------------------------------------------------------------------

% Construct a cell array storing the names of model parameters.
ParamsNames = {'P0(0)','P1(0)','P2(0)','Lambda1','Lambda2','Theta1','Theta2','Delta','Gamma','K'};
% Get the number of model parameters that are allowed to vary.
varying_parameters_num = length(ParamIndices);
% Get the optimal values for T1 and T2.
T1 = Topts(:,1);
T2 = Topts(:,2);
% Get the optimal values for S1 and S2.
S1 = Sopts(:,2);
S2 = Sopts(:,3);
if(varying_parameters_num==2)
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
    % Store separately the values and names of the parameters that are
    % allowed to vary.
    VaryingParamIndex = ParamIndices(1);
    ThetaParamIndex = ParamIndices(2);
    VaryingParam = Params(:,VaryingParamIndex);
    VaryingParamName = ParamsNames{VaryingParamIndex};
    ThetaParam = Params(:,ThetaParamIndex);
    ThetaParamName = ParamsNames{ThetaParamIndex};
    % Retrieve the actual range of the varying parameter.
    VaryingParamRange = unique(VaryingParam);
    % Initialize vector storing the optimal influence ranges.
    OptimalInfluenceRanges = [];
    % For each distinct value of the varying parameter collect the values
    % of T1_opt,T2_opt,S1_opt,S2_opt,Flags1 and Flags2 in order to
    % determine the influence regions of interest that are defined below:
    % Ra(Theta) = {Theta: T1_opt(Theta)>=T2_opt(Theta) and S1_opt(Theta)<=S2_opt(Theta)}
    % Rb(Theta) = {Theta: T2_opt(Theta)>=T1_opt(Theta) and S2_opt(Theta)<=S1_opt(Theta)}
    for varying_param_value = VaryingParamRange'
        varying_param_indices = find(VaryingParam==varying_param_value);
        t1_opt = T1(varying_param_indices);
        t2_opt = T2(varying_param_indices);
        s1_opt = S1(varying_param_indices);
        s2_opt = S2(varying_param_indices);
        %Ra_indices = find(and((t1_opt>t2_opt),s1_opt<s2_opt));
        Ra_indices = find(and((t1_opt-t2_opt>0.001),s1_opt-s2_opt<0.001));
        Ra_theta_indices = varying_param_indices(Ra_indices);
        Ra_thetas = ThetaParam(Ra_theta_indices);
        %Rb_indices = find(and((t2_opt>t1_opt),s2_opt<s1_opt));
        Rb_indices = find(and((t2_opt-t1_opt>0.001),s2_opt-s1_opt<0.001));
        Rb_theta_indices = varying_param_indices(Rb_indices);
        Rb_thetas = ThetaParam(Rb_theta_indices);
        R_thetas = union(Ra_thetas,Rb_thetas);
        if(~isempty(R_thetas))
            optimal_influence_range = max(R_thetas) - min(R_thetas);
        else
            optimal_influence_range = 0;
        end;
        OptimalInfluenceRanges = [OptimalInfluenceRanges;optimal_influence_range];
    end;
    % Generate the title string for the figure.
    TitleString = strcat([Const1Name ' = ' Const1ValueString ' ' Const2Name ' = ' ...
                        Const2ValueString ' ' Const3Name ' = ' Const3ValueString ' ' ...
                        Const4Name ' = ' Const4ValueString '\n' Const5Name ' = ' Const5ValueString ' '...
                        Const6Name ' = ' Const6ValueString ' ' Const7Name ' = ' Const7ValueString ' ' ...
                        Const8Name ' = ' Const8ValueString]);
    TitleName = sprintf(TitleString);
    % Construct figure.
    FigureName = strcat(['Influence Range  with respect to ' VaryingParamName]);
    figure('Name',FigureName);
    hold on
    plot(VaryingParamRange,OptimalInfluenceRanges,'*k','LineWidth',1.4);
    xlabel(VaryingParamName);
    ylabel('Less Influential Investment Region');
    grid on
    title(TitleName);
end;




end

