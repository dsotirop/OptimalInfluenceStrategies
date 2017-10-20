function PlotBestResponseCurves(Ra,Rb)

% This function plots the best response curves stored as point pairs within
% matrices Ra and Rb.
% Ra: stores point pairs of the form (TAopt,TB) for any given TB.
% Rb: stotes point pairs of the form (TA,TBopt) for any given TA.

% Given that the values of TA and TB as independent variables are
% restricted with the [0,1] interval the actual limits will be obtain by
% measuring the boundary values for TAopt and TBopt.

% Initially obtain the individual vectors TA, TB, TAopt and TBopt.
TA = Rb(:,1);
TB = Ra(:,2);
TAopt = Ra(:,1);
TBopt = Rb(:,2);

% Get the actual limits for TAopt and TBopt.
TAmin = min(TAopt);
TAmax = max(TAopt);
TBmin = min(TBopt);
TBmax = max(TBopt);

% Set the number of possible equilibrium points.
Neq = 6;

% Compute all the pairwise Manhattan distances between points in Ra and Rb.
Dab = mandist(Ra,Rb');

% Get the number of pairs of points between the two curves.
Nab = numel(Dab);

% Derive the vector of sorted distances in ascending order.
Dab_sort = sort(reshape(Dab,1,Nab),'ascend');

% Get the miminum and maximum distance values within Dab.
Dab_min = Dab_sort(1);
Dab_max = Dab_sort(end);

% Get the top Neq minimum distances for pairs of points between Ra and Rb
% which will serve as candidates for equilibrium points.
Deq = Dab_sort(1:Neq);


% Create new plotting figure for the pairwise distances.
% figure('Name','Pairwise Distances');
% clf
% plot([1:1:Nab],Dab_sort,'.b','LineWidth',2);
% grid on

% Initialize indices vectors for storing the positions in Dab that
% correspond to the distances in Deq.
IJeq = [];
Peq = [];
for k = 1:1:Neq
    fprintf('Number %d closest distance: %d\n',k,Deq(k));
    [i,j] = find(Dab==Deq(k));
    IJeq = [IJeq;[i j]];
    Peq = [Peq;[Ra(i,1),Ra(i,2)]];
    Peq = [Peq;[Rb(j,1),Rb(j,2)]];
end;

% Create new plotting figure for the best response curves.
figure('Name','Simplified Oligopolistic Game Best Responses');
clf
hold on
% Plot the individual reaction functions.
plot(TAopt,TB,'.r','LineWidth',2);
plot(TA,TBopt,'.g','LineWidth',2);
% Plot the candidate equilibrium points.
plot(Peq(:,1),Peq(:,2),'*k','LineWidth',1);
% Set axis limits.
axis([TAmin max(TAmax,1) TBmin max(TBmax,1)]);
grid on
xlabel('T_A Optimal');
ylabel('T_B Optimal');
% Shade the area of accepted solutions.
% alpha(0.1);
% X = [0 1];
% Y = [1 0];
% h = area(X,Y);
% h.FaceColor = (1/255)*[176,224,230];
end

