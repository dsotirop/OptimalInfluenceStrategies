% -------------------------------------------------------------------------
% Define a function for plotting the triplets of model variables with
% respect to the varying parameter.
function PlotTripletsX(Triplet,BaseComparison)
% -------------------------------------------------------------------------
% Triplet: is a structure storing all the necessary information for carrying
%          out required plots.
% BaseComparison: is a boolean variable that controls whether the minus and
%                 plus experiments will be plotted along with the base 
%                 experiment for comparison purposes.    
% -------------------------------------------------------------------------
% Define the rgb colors to be used for plotting the variables of each
% experiment.
RedColor   =  [(255/255),(0/255),(0/255)];
GreenColor =  [(0/255),(128/255),(0/255)];
BlueColor  =  [(0/255),(0/255),(255/255)];

% Set color for firm A plots.
ColorA = RedColor;
ColorB = GreenColor;
ColorC = BlueColor;
% -------------------------------------------------------------------------
% Loop through the various model variables that need to be plotted.
for variable_index = 1:1
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
    Lo = length(Triplet.TitleConstParams);
    for t = 1:Lo
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
% -------------------------------------------------------------------------   
    % Plot the optimal variables for the two firms with respect to the
    % varying parameter. Mind that when the optimal limiting influences are
    % plot an additional plot for the limiting influence of the consumer
    % should be incorporated in the graph.
    % Set the independent variables of each plot appearing on the x-axis.
    x_minus = Triplet.MinusVaryingParam;
    x_base = Triplet.BaseVaryingParam;
    x_plus = Triplet.PlusVaryingParam;
    % Set the dependent variables of each plot appearing on the y-axis.
    y_minus = Triplet.MinusExperimentVariables{variable_index};
    y_base = Triplet.BaseExperimentVariables{variable_index};
    y_plus = Triplet.PlusExperimentVariables{variable_index};
% -------------------------------------------------------------------------    
    % Do the actual plotting for the minus experiment for the two firms.
    L = length(Triplet.LegendConstParamsNames);
    % Set the third line of the title string for the current figure.
    for t = 1:L
        if(t==1)
            ThirdLineTitleString = sprintf('%s = %s',...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.MinusLegendConstParams(t)));
        else
            ThirdLineTitleString = sprintf('%s  %s = %s',...
                                    ThirdLineTitleString,...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.MinusLegendConstParams(t)));
        end
    end
    subplot(3,1,1);
    hold on
    if(variable_index ~= 2)
        p_minus_a = plot(x_minus,y_minus(:,1),'-','LineWidth',2.0);
        p_minus_b = plot(x_minus,y_minus(:,2),'-','LineWidth',2.0);
        set(p_minus_a,'Color',ColorA);
        set(p_minus_b,'Color',ColorB);
    else
        p_minus_a = plot(x_minus,y_minus(:,1),'-','LineWidth',2.0);
        p_minus_b = plot(x_minus,y_minus(:,3),'-','LineWidth',2.0);
        set(p_minus_a,'Color',ColorA);
        set(p_minus_b,'Color',ColorB);
        p_minus_c = plot(x_minus,y_minus(:,2),'-','LineWidth',2.0);
        set(p_minus_c,'Color',ColorC);
    end
    % -----------------------------------------------------------------
    % Add the option to present the base experiment in the same graph.
    if(BaseComparison)
        if(variable_index ~= 2)
            p_base_a = plot(x_base,y_base(:,1),':','LineWidth',2.0);
            p_base_b = plot(x_base,y_base(:,2),':','LineWidth',2.0);
            set(p_base_a,'Color',ColorA);
            set(p_base_b,'Color',ColorB);
        else
            p_base_a = plot(x_base,y_base(:,1),':','LineWidth',2.0);
            p_base_b = plot(x_base,y_base(:,3),':','LineWidth',2.0);
            set(p_base_a,'Color',ColorA);
            set(p_base_b,'Color',ColorB);
            p_base_c = plot(x_base,y_base(:,2),':','LineWidth',2.0);
            set(p_base_c,'Color',ColorC);
        end
    % Adjust the axes limits.
    x_min = min([min(x_minus) min(x_base)]);
    x_max = max([max(x_minus) max(x_base)]);
    y_min = min([min(min(y_minus)) min(min(y_base))]);
    y_max = max([max(max(y_minus)) max(max(y_base))]);
    if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
        y_min = 0;
        y_max = 1;
    end
    axis([x_min x_max y_min y_max]);
    end    
    % -----------------------------------------------------------------
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString,SecondLineTitleString,...
                   ThirdLineTitleString}; 
    % Set the title for the current figure.
    title(TitleString,'fontsize',12);
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',12);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',12);
    % Enable plotting grid.
    grid on
    % Set the axes limits in case the BaseComparison is set to false.
    if(~BaseComparison)
        x_min = min(x_plus);
        x_max = max(x_plus);
        y_min = min(min(y_plus));
        y_max = max(max(y_plus));
        if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
            y_min = 0;
            y_max = 1;
        end
        axis([x_min x_max y_min y_max]);
    end
    hold off
% -------------------------------------------------------------------------    
    % Do the actual plotting for the base experiment for the two firms.
    for t = 1:L
        if(t==1)
            ThirdLineTitleString = sprintf('%s = %s',...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.BaseLegendConstParams(t)));
        else
            ThirdLineTitleString = sprintf('%s  %s = %s',...
                                   ThirdLineTitleString,...
                                   Triplet.LegendConstParamsNames{t},...
                                   num2str(Triplet.BaseLegendConstParams(t)));
        end
    end
    subplot(3,1,2);
    hold on
    if(variable_index ~= 2)
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.0);
        p_base_b = plot(x_base,y_base(:,2),'-','LineWidth',2.0);
        set(p_base_a,'Color',ColorA);
        set(p_base_b,'Color',ColorB);
    else
        p_base_a = plot(x_base,y_base(:,1),'-','LineWidth',2.0);
        p_base_b = plot(x_base,y_base(:,3),'-','LineWidth',2.0);
        set(p_base_a,'Color',ColorA);
        set(p_base_b,'Color',ColorB);
        p_base_c = plot(x_base,y_base(:,2),'-','LineWidth',2.0);
        set(p_base_c,'Color',ColorC);
    end
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString,SecondLineTitleString,...
                   ThirdLineTitleString};
    % Set the title for the current figure.
    title(TitleString,'fontsize',12);
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',12);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',12);
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
    hold off
% -------------------------------------------------------------------------    
    % Do the actual plotting for the plus experiment for the two firms.
    for t = 1:L
        if(t==1)
            ThirdLineTitleString = sprintf('%s = %s',...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.PlusLegendConstParams(t)));
        else
            ThirdLineTitleString = sprintf('%s  %s = %s',...
                                    ThirdLineTitleString,...
                                    Triplet.LegendConstParamsNames{t},...
                                    num2str(Triplet.PlusLegendConstParams(t)));
        end
    end
    subplot(3,1,3);
    hold on
    if(variable_index ~= 2)
        p_plus_a = plot(x_plus,y_plus(:,1),'-','LineWidth',2.0);
        p_plus_b = plot(x_plus,y_plus(:,2),'-','LineWidth',2.0);
        set(p_plus_a,'Color',ColorA);
        set(p_plus_b,'Color',ColorB);
    else
        p_plus_a = plot(x_plus,y_plus(:,1),'-','LineWidth',2.0);
        p_plus_b = plot(x_plus,y_plus(:,3),'-','LineWidth',2.0);
        set(p_plus_a,'Color',ColorA);
        set(p_plus_b,'Color',ColorB);
        p_plus_c = plot(x_plus,y_plus(:,2),'-','LineWidth',2.0);
        set(p_plus_c,'Color',ColorC);
    end
    % -----------------------------------------------------------------
    % Add the option to present the base experiment in the same graph.
    if(BaseComparison)
        if(variable_index ~= 2)
            p_base_a = plot(x_base,y_base(:,1),':','LineWidth',2.0);
            p_base_b = plot(x_base,y_base(:,2),':','LineWidth',2.0);
            set(p_base_a,'Color',ColorA);
            set(p_base_b,'Color',ColorB);
        else
            p_base_a = plot(x_base,y_base(:,1),':','LineWidth',2.0);
            p_base_b = plot(x_base,y_base(:,3),':','LineWidth',2.0);
            set(p_base_a,'Color',ColorA);
            set(p_base_b,'Color',ColorB);
            p_base_c = plot(x_base,y_base(:,2),':','LineWidth',2.0);
            set(p_base_c,'Color',ColorC);
        end
    % Adjust the axes limits.
    x_min = min([min(x_plus) min(x_base)]);
    x_max = max([max(x_plus) max(x_base)]);
    y_min = min([min(min(y_plus)) min(min(y_base))]);
    y_max = max([max(max(y_plus)) max(max(y_base))]);
    if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
        y_min = 0;
        y_max = 1;
    end
    axis([x_min x_max y_min y_max]);    
    end    
    % -----------------------------------------------------------------
    % Set the title string for the current figure.
    TitleString = {FirstLineTitleString,SecondLineTitleString,...
                   ThirdLineTitleString};
    % Set the title for the current figure.
    title(TitleString,'fontsize',12);
    % Set the x-label for the current figure.
    xlabel(XLabelString,'fontsize',12);
    % Set the y-lable for the current figure.
    ylabel(YLabelString,'fontsize',12);
    % Enable plotting grid.
    grid on
    % Set the axes limits in case the BaseComparison is set to false.
    if(~BaseComparison)
        x_min = min(x_plus);
        x_max = max(x_plus);
        y_min = min(min(y_plus));
        y_max = max(max(y_plus));
        if(y_max-y_min<10^-Triplet.MinimumDigitsAccuracy)
            y_min = 0;
            y_max = 1;
        end
        axis([x_min x_max y_min y_max]);
    end
    hold off
    % Set the position and the size of the figure.
    set(f,'units','points','position',[10,10,600,600])
end
% -------------------------------------------------------------------------
end