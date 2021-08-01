% This scipt file performs all necessary derivations with respect to the
% limiting influences Sj and S_j that appear in the Simplified
% Oligopolistic Influence Model.

clc
clear

% Define fundamental symbolic quantities.
syms Tj T_j Lj L_j
syms Sj(Tj,T_j) S_j(Tj,T_j)

% Define the quantities for limiting influences.
Sj = L_j * Tj / (Lj*L_j + L_j*Tj + Lj*T_j);
S_j = Lj * T_j / (Lj*L_j + L_j*Tj + Lj*T_j);

% Compute the first derivatives of the quantity Sj with respect to Tj and
% T_j.
DSjTj = diff(Sj,Tj);
DSjTj = simplify(collect(expand(DSjTj)));
DSjT_j = diff(Sj,T_j);
DSjT_j = simplify(collect(expand(DSjT_j)));
% Compute the second derivatives of the quantities DSjTj and DSjT_j with
% respect to Tj and T_j.
D2SjTjTj = diff(DSjTj,Tj);
D2SjTjTj = simplify(collect(expand(D2SjTjTj)));
D2SjT_jT_j = diff(DSjT_j,T_j);
D2SjT_jT_j = simplify(collect(expand(D2SjT_jT_j)));
D2SjTjT_j = diff(DSjTj,T_j);
D2SjTjT_j = simplify(collect(expand(D2SjTjT_j)));

% Compute the first derivatives of the quantity S_j with respect to Tj and
% T_j.
DS_jTj = diff(S_j,Tj);
DS_jTj = simplify(collect(expand(DS_jTj)));
DS_jT_j = diff(S_j,T_j);
DS_jT_j = simplify(collect(expand(DS_jT_j)));
% Compute the second derivatives of the quantities DS_jTj and DS_jT_j with
% respect to Tj and T_j.
D2S_jTjTj = diff(DS_jTj,Tj);
D2S_jTjTj = simplify(collect(expand(D2S_jTjTj)));
D2S_jT_jT_j = diff(DS_jT_j,T_j);
D2S_jT_jT_j = simplify(collect(expand(D2S_jT_jT_j)));
D2S_jTjT_j = diff(DS_jTj,T_j);
D2S_jTjT_j = simplify(collect(expand(D2S_jTjT_j)));