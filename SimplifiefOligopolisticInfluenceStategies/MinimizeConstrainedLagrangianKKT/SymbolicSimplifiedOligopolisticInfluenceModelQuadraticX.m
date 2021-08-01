% This script file performs fundamental symbolic computations in order to
% derive the equations governing a simplified oligopolistic influence
% strategies model concering a network G = {F1,F2,...,Fn,C} of (n) firms 
% F = {F1,F2,...,Fn} and one (1) consumer C.

% Clear command window and workspace.
clc
clear

% Set the number n of firms pertaining to the network.
n = 3;

% -------------------------------------------------------------------------
% Setup the network-related symbolic variables:
% -------------------------------------------------------------------------
% Create the symbolic variables capture the direct influnce exerted 
% by consumer C on each firm in F = {F1,F2,...,Fn}.
L = sym('L%d',[1 n]);
% Create the symbolic variables correspond to the fundamental optimization 
% variables of the model capturing the levels of direct influence that each 
% firm exerts on the consumer C.
T = sym('T%d',[1 n]);