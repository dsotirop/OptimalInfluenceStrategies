function PlotFundamentalFunctions3D(C,G,LA,LB,PA,PB,alpha,beta,gamma,TRange)

% This function plots the following optimization model quantities
% (i):   Firm A Profit
% (ii):  Firm B Profit
% (iii): Firm A Profit First Derivative
% (iv):  Firm B Profit First Derivative
% (v):   Firm A Profit Second Derivative
% (vi):  Firm B Profit Second Derivative
% on an meshgrid which is synthesized of (TA,TB) value pairs such that:
% (TA,TB) in TRange x TRange. TRange is a vector of the form [tmin:dt:tmax]
% which controls the range values for each individual optimization variable.

% Construct the meshgrid of (TA,TB) value pairs for the given range vector.
[TA,TB] = meshgrid(TRange,TRange);

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
set(ha,'FaceColor',[1 0 0], 'FaceAlpha',0.9,'EdgeAlpha', 0);
set(hb,'FaceColor',[0 1 0], 'FaceAlpha',0.9,'EdgeAlpha', 0);
xlabel('T_A');
ylabel('T_B');
zlabel('f_a and f_b');
legend({'f_a','f_b'});

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
zlabel('D^2_a');
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


end

