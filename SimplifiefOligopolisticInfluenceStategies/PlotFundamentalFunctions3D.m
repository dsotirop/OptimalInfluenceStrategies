function PlotFundamentalFunctions3D(C,G,LA,LB,PA,PB,alpha,beta,gamma,TRange)

% This function plots the following optimization model quantities
% (i):    Firm A Profit
% (ii):   Firm B Profit
% (iii):  Firm A Profit First Derivative
% (iv):   Firm B Profit First Derivative
% (v):    Firm A Profit Second Derivative
% (vi):   Firm B Profit Second Derivative
% (vii):  Firm A Revenue
% (viii): Firm B Revenue
% (ix):   Firm A Revenue First Derivative
% (x):    Firm B Revenue First Derivative
% (xi):   Firm A Revenue Second Derivative
% (xii):  Firm B Revenue Secondn Derivative
% on an meshgrid which is synthesized of (TA,TB) value pairs such that:
% (TA,TB) in TRange x TRange. TRange is a vector of the form [tmin:dt:tmax]
% which controls the range values for each individual optimization variable.

% Define the grayscale color palette.
BlackColor = [0,0,0];
LightsLateGrayColor = [(119/255),(136/255),(153/255)];
GrayColor = [(128/255),(128/255),(128/255)];
GainsboroColor = [(220/255),(220/255),(220/255)];
DarksLateGrayColor = [(47/255),(79/255),(79/255)];
CrimsonColor = [(220/255),(20/255),(60/255)];
DeepSkyBlueColor = [(0/255),(191/255),(255/255)];
LimeColor = [(0/255),(255/255),(0/255)];
DimGrayColor = [(105/255),(105/255),(105/255)];
BlueVioletColor = [(138/255),(43/255),(226/255)];

% Construct the meshgrid of (TA,TB) value pairs for the given range vector.
[TA,TB] = meshgrid(TRange,TRange);

% -------------------------------------------------------------------------
%                PROFIT FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the profit functions for
% the two firms evaluated at each point of the previously defined mesh
% grid. In particular, profit functions will be plotted seperately and
% combined in a single plot.
fa = FirmAProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
fb = FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Profit');
mesh(TA,TB,fa);
xlabel('T_A');
ylabel('T_B');
zlabel('f_a');
figure('Name','Firm B Profit');
mesh(TA,TB,fb);
xlabel('T_A');
ylabel('T_B');
zlabel('f_b');
figure('Name','Firms A and B Profits');
ha = surf(TA,TB,fa);
hold on
hb = surf(TA,TB,fb);
% Set the red color for the surface of fa and the green color for the 
% surface of fb. 
set(ha,'FaceColor',BlueVioletColor,'FaceAlpha',0.9,'EdgeAlpha', 0);
set(hb,'FaceColor',CrimsonColor,'FaceAlpha',0.9,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('f_a and f_b');
legend({'f_a','f_b'});
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%                REVENUE FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the revenue functions for
% the two firms evaluated at each point of the previously defined mesh
% grid. In particular, profit functions will be plotted seperately and
% combined in a single plot.
Ra = FirmARevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Rb = FirmBRevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Revenue');
mesh(TA,TB,Ra);
xlabel('T_A');
ylabel('T_B');
zlabel('R_a');
figure('Name','Firm B Revenue');
mesh(TA,TB,Rb);
xlabel('T_A');
ylabel('T_B');
zlabel('R_b');
figure('Name','Firms A and B Revenues');
ha = surf(TA,TB,Ra);
hold on
hb = surf(TA,TB,Rb);
% Set the red color for the surface of ra and the green color for the 
% surface of rb. 
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.9,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.9,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('R_a and R_b');
legend({'R_a','R_b'});
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%           FIRST DERIVATIVE PROFIT FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the first derivative profit 
% functions for the two firms evaluated at each point of the previously 
% defined mesh grid. In particular, profit functions will be plotted 
% seperately and combined in a single plot.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Profit First Derivative');
mesh(TA,TB,Da);
xlabel('T_A');
ylabel('T_B');
zlabel('D_a');
figure('Name','Firm B Profit First Derivative');
mesh(TA,TB,Db);
xlabel('T_A');
ylabel('T_B');
zlabel('D_b');
figure('Name','Firms A and B Profits First Derivatives');
ha = surf(TA,TB,Da);
hold on
hb = surf(TA,TB,Db);
% Set the red color for the surface of Da and the green color for the 
% surface of Db. 
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('D_a and D_b');
legend({'D_a','D_b'});
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%           FIRST DERIVATIVE REVENUE FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the first derivative revenue 
% functions for the two firms evaluated at each point of the previously 
% defined mesh grid. In particular, profit functions will be plotted 
% seperately and combined in a single plot.
DRa = FirmARevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
DRb = FirmBRevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Revenue First Derivative');
mesh(TA,TB,DRa);
xlabel('T_A');
ylabel('T_B');
zlabel('DR_a');
figure('Name','Firm B Revenue First Derivative');
mesh(TA,TB,DRb);
xlabel('T_A');
ylabel('T_B');
zlabel('DR_b');
figure('Name','Firms A and B Revenues First Derivatives');
ha = surf(TA,TB,DRa);
hold on
hb = surf(TA,TB,DRb);
% Set the red color for the surface of DRa and the green color for the 
% surface of DRb. 
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('DR_a and DR_b');
legend({'DR_a','DR_b'});
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%           SECOND DERIVATIVE PROFIT FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the second derivative profit 
% functions for the two firms evaluated at each point of the previously 
% defined mesh grid. In particular, profit functions will be plotted 
% seperately and combined in a single plot.
DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Profit Second Derivative');
mesh(TA,TB,DDa);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2_a');
figure('Name','Firm B Profit Second Derivative');
mesh(TA,TB,DDb);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2_b');
figure('Name','Firms A and B Profits Second Derivatives');
ha = surf(TA,TB,DDa);
hold on
hb = surf(TA,TB,DDb);
% Set the red color for the surface of DDa and the green color for the 
% surface of DDb. 
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2_a and D^2_b');
legend({'D_a','D_b'});
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%           SECOND DERIVATIVE REVENUE FUNCTIONS 3 DIMENSIONAL PLOTS:
% -------------------------------------------------------------------------
% Construct the 3 dimensional surface plots of the second derivative profit 
% functions for the two firms evaluated at each point of the previously 
% defined mesh grid. In particular, profit functions will be plotted 
% seperately and combined in a single plot.
DDRa = FirmARevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
DDRb = FirmBRevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
figure('Name','Firm A Revenue Second Derivative');
mesh(TA,TB,DDRa);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2R_a');
figure('Name','Firm B Revenue Second Derivative');
mesh(TA,TB,DDRb);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2R_b');
figure('Name','Firms A and B Revenue Second Derivatives');
ha = surf(TA,TB,DDRa);
hold on
hb = surf(TA,TB,DDRb);
% Set the red color for the surface of DDRa and the green color for the 
% surface of DDRb. 
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.8,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('D^2R_a and D^2R_b');
legend({'D^2R_a','D^2R_b'});
% -------------------------------------------------------------------------

end

