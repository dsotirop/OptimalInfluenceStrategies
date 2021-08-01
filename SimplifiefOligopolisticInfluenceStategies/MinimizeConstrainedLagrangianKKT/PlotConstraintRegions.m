function PlotConstraintRegions(alpha,beta,LA,LB,PA,PB,RSOL)

% This function plots the constraint regions for the non-linear
% optimization problem that underlies the continuous two player game for
% the Simplified Oligopolistic Optimal Influence Model. Mind that the
% aforementioned constraint regions are plotted only for the case where the
% competition related parameter beta is negative which activates a set of
% 3 linear constraints in order to ensure the positivity of the optimal
% demand functions. Moreover, it should be noted that the set of linear
% constraints is obtained for the special case where the marginal cost
% parameter c is set to zero (c = 0).

% RSOL stores the vector of optimal solutions for the investment levels for
% the two firms in the case where the non-linear optimization process has 
% reached an optimal solution that satisfies all the necessary constraints.  

% Consider the three linear equations that define the corresponding constraints
% in the form g_j(TA,TB) = a_j * TA + b_j * TB + c_j = 0.
% Re-express the previous equations in the form TB_j = f_j(TA) so that the 
% three lines may be plotted as a sequence of points (TA(k),TB_j(k)). 
TA = 0:0.01:1;
% Set the 1st linear constraint equation: g1(TA,TB) = 0.
TB_1 = 1 - TA; 
% Set the 2nd linear constraint equation: g2(TA,TB) = 0.
TB_2 = (-alpha*LB*TA-LA*LB*(alpha*PA+beta*PB)) / (beta*LA);
% Set the 3rd linear constraint equation: g3(TA,TB) = 0.
TB_3 = (-beta*LA*TA-LA*LB*(alpha*PB+beta*PA)) / (alpha*LB);

% Initialize the new plotting window.
figure('Name','Linear Constraint Regions');
%axis([0 1 0 1]);
hold on
% Plot the main axis of the [0,1] x [0,1] interval.
plot(TA,zeros(1,length(TA)),'k','LineWidth',2.0);
plot(zeros(1,length(TA)),TA,'k','LineWidth',2.0);
% Plot the 1st linear equation.
plot(TA,TB_1,'-k','LineWidth',2.0);
% Plot the 2nd linear equation.
plot(TA,TB_2,'-r','LineWidth',2.0);
% Plot the 3rd linear equation.
plot(TA,TB_3,'-b','LineWidth',2.0);

% Plot the area satisfying the 1st linear consraint equation: g1(TA,TB) < 0.
P1 = patch([TA fliplr(TA)],[TB_1 min(ylim)*ones(size(TB_1))],'k');
P1.FaceAlpha = 0.3;
% Plot the area satisfying the 2nd linear constraint equation: g2(TA,TB) > 0.
P2 = patch([TA fliplr(TA)],[TB_2 min(ylim)*ones(size(TB_2))],'r');
P2.FaceAlpha = 0.3;
% Plot the area satisfying the 3rd linear constraint equation: g3(TA,TB) > 0.
P3 = patch([TA fliplr(TA)],[TB_3 max(ylim)*ones(size(TB_3))],'b');
P3.FaceAlpha = 0.3;
% Plot the solution point if it exists.
if(~isempty(RSOL))
    plot(RSOL(1),RSOL(2),'*m','LineWidth',2.0);
end
hold off
grid on

end