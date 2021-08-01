% This script file facilitates the conduction of the fundamental symbolic
% operations that underly the modified version of the Simplified 
% Oligopolistic Influence Model. The main diversification of this approach
% lies upon the utilization of an extended version of the 
% "Quality and Competition" model which was originally proposed by Banker 
% et al.

% For this particular formulation we utilize the original version of the
% Demand function given by:
%
% Qx = (Xx - K*Xy) - (Px - M*Py) [1]  
% Qy = (Xy - K*Xx) - (Py - M*Px) [2]
% 
% where Xx and Xy are the perceived qualities of the products produced by
% firms x and y. Moreover, quantities Xx and Xy correspond to the consumer's
% limiting beliefs about the two products. Px and Py are the price  levels  
% chosen by firms x and y, whereas, M is the degree of competition on the 
% and K is the degree of competition on the quality level.

% Letting Yx and Yy be the true qualitues of the products produced by firms
% x and y, the evenue functions for the two firms will be given as:
%
%    Rx = (Px - R*Yx) * Qx [3]
%    Ry = (Py - R*Yy) * Qy [4]
%
% where R is the quality-related cost parameter.

% Cost functions Cx and Cy which will be defined according to Banker et al.
% as:
%               2                    
%    Cx = G * Xx   [5]
%               2 
%    Cy = G * Xy   [6]

% Clear command window and workspace.
clc
clear

% Define the required symbolic variables.
syms Xx Xy Px Py Yx Yy G
% Define additional symbolic variables.
syms px py Sx Sy So
% Define additional symbolic variables.
syms Lx Ly Tx Ty

% Define the expressions for the quantities.
Qx = (Xx - K*Xy) - (Px - M*Py);
Qy = (Xy - K*Xx) - (Py - M*Px);
% Define the expressions for the revenues.
Rx = (Px - R*Yx) * Qx;
Ry = (Py - R*Yy) * Qy;
% Define the expressions for the cost functions.
Cx = G * Xx^2;
Cy = G * Xy^2;
% Define the expressions for the profits.
Fx = Rx - Cx;
Fy = Ry - Cy;