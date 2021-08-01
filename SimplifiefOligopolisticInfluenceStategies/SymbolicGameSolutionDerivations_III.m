% This script file investigates the possibilty of acquring an analytical
% solution for the continuous game that underlies the Simplified
% Oligopolistic Optimal Influence model.

% This analysis reles on the following derivations for the profit functions
% Fx and Fy of the two firms A and B that interact with consumer C through 
% the network:
%
%        2         2
% Fx = Qx  - G * Tx   [1] (Profix function for Firm x) 
%        2         2 
% Fy = Qy -  G * Ty   [2] (Profix function for Firm y)

% Qx and Qy are the second stage equilibrium demand functions for firms A
% and B that may be expressed with respect to the influence levels Sx and 
% Sy as:
%
% Qx = Axy * Sx + Byx * Sy + Gxy [3]
% Qy = Ayx * Sy + Bxy * Sx + Gyx [4]

% The influnce levels Sx and Sy of the two firms in the network are given
% as:
%       Ly * Tx                Lx * Ty
% Sx = --------- [5] and Sy = --------- [6] where
%         So                     So
%
% So = Lx * Ly + Ly * Tx + Lx * Ty [7]

% Equations [5] and [6] may be utilized in order to express the direct
% influence levels of the two firms as functions of the corresponding
% influce levels:
%
%         Lx * Sx                   Ly * Sy   
% Tx = ------------- [8] and Ty = ------------- [9]
%       1 - Sx - Sy                1 - Sx - Sy  
%
% In this context, equations [8] and [9] may be used to express So in the
% following way:
%                      Lx * Ly
%               So = ----------- [10]
%                    1 - Sx - Sy

% The First Order Conditions for the underlying optimization problems imply
% that:
%
%           dFx              dQx
%  DFxTx = ----- = 2 * Qx * ----- - 2 * G * Tx == 0 [11]
%           dTx              dTx
%
%           dFy              dQy
%  DFyTx = ----- = 2 * Qy * ----- - 2 * G * Ty == 0 [12]
%           dTy              dTy

% The first order derivatives of the equilibrium demand functions Qx and Qy
% with respect to Tx and Ty respectively may be expressed in the following
% way:
%
%           dQx       Ly
% DQxTx = -------  = ---- * (a - Qx) [13]
%           dTx       So
%
%           dQy       Lx
% DQyTy = -------  = ---- * (a - Qy) [13]
%           dTy       So

% Clear command window and workspace.
clc
clear 

% Define the required symbolic variables.
syms Lx Ly Tx Ty Sx Sy So a b Axy Ayx Bxy Byx Gxy Gyx G
% Define the expressions for the equilibrium demand functions.
Qx = Axy * Sx + Byx * Sy + Gxy;
Qy = Ayx * Sy + Bxy * Sx + Gyx;
% Define the expressions for the direct influence levels.
Tx = (Lx * Sx)/(1 - Sx - Sy);
Ty = (Ly * Sy)/(1 - Sx - Sy);
% Define the expression for the quantity So.
So = (Lx * Ly)/(1 - Sx - Sy);
% Define the expressions for the first order derivatives of Qx and Qy with
% respect to Tx and Ty respectively.
DQxTx = (Ly/So) * (a - Qx);
DQyTy = (Lx/So) * (a - Qy);
% Define the expressions for the system of equations that define the Nash
% equilibrium of the game.
DFxTx = Qx * DQxTx - G * Tx;
DFyTy = Qy * DQyTy - G * Ty;
% Express quantities DFxTx and DFyTy as fractions.
[DFxTx_n,DFxTx_d] = numden(DFxTx);
[DFyTy_n,DFyTy_d] = numden(DFyTy);
% Simplify the expressions for the numerators.
DFxTx_n = simplify(collect(expand(DFxTx_n)));
DFyTy_n = simplify(collect(expand(DFyTy_n)));
% Define the system of equations that need to be solved.
% Sys = [DFxTx_n==0;DFyTy_n==0];
% Sol = solve(Sys,[Sx Sy])