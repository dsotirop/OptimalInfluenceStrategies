% This script file performs fundamental symbolic computations in order to
% derive the equations governing a simplified oligopolistic influence
% strategies model concering a network G = {A,C,B} of two firms {A,B} and
% one consumer C.

% This script file is an extension to the original 
% "SymbolicSimplifiedInfluenceModel.m" file with the only differentation
% being the addition of a quadratic cost term on the expressions for the
% profits of the two firms A and B. Thus, the expression for the two
% profits function will be subsequently defined as:
% (i):  FA_opt = (pA_opt - C)^2 - G * TA^2
% (ii): FB_opt = (pB_opt - C)^2 - G * TB^2

% Clear command window and workspace.
clc
clear

% Setup the network-related symbolic variables.
syms LA LB % These symbolic variables capture the direct influnce exerted 
           % by consumer C on the two firms A and B.
syms TA TB % These symbolic variables correspond to the fundamental 
           % optimization variables of the model capturing the levels of 
           % direct influence that each firm exerts on the consumer C.
syms VSA VSC VSB % These are auxiliary symbolic varibles which are utilized 
                 % in order to compute the limiting influence vector
                 % S = [SA SC SB] when the interaction amonst the network
                 % agents will reach a consensial belief XA and XB for each
                 % product A and B.

% The naming convention for the symbolic variables holding the intial
% agents' beliefs is the following:
% PXY is the initial belief of firm X in {A,B} for product Y in {A,B}.
% PZ is the initial belief of consumer X for product Z in {A,B}.
syms PAA PA PBA % These symbolic variables correspond to the initial beliefs
                % the network agents hold for the product A such that:
                % (i): PAA the initial belief firm A holds for product A.
                % (ii): PA the initial belief consumer C holds for product A.
                % (iii): PBA the initial belief firm B holds for product A.
syms PAB PB PBB % These symbolic variables correspond to the initial beliefs
                % the network agents hold for the product B such that:
                % (i): PAB the initial belief firm A holds for product B.
                % (ii): PB the initial belief consumer C holds for product B.
                % (iii): PBB the initial belief firm B holds for product B.

syms XA XB % These symbolic variables correspond to the limiting beliefs of 
           % each agent in the network for products A and B.
syms QA QB % These symbolic variables correspond to the quantities consumed 
           % by consumer C from each product A and B.
syms pA pB % These symbolic variables correspond to the prices for each 
           % product A and B.
syms M K % These symbolic variables correspond to the sensitivity coefficients.
syms C % This symbolic variable corresponds to the marginal cost.
syms G % This symbolic variable corresponds to the marginal influence cost or Gamma.

% Setup social interaction matrix components.
T = [1-LA LA 0;TA 1-TA-TB TB;0 LB 1-LB];
% Setup the components of the final influence vector.
VS = [VSA VSC VSB];
% Setup the corresponding system initiating from the left eigenvalues 
% computation problem. Mind that the number of equations required from Y==0
% is 2. In this case we take into account equations originating from rows
% (1) and (3) as indicated by our analytical derivations.
Y = VS * (T - eye(3));
System = [Y(1)==0;Y(3)==0;VSA+VSC+VSB==1];
Solution = solve(System,VSA,VSC,VSB);
SA = Solution.VSA;
SC = Solution.VSC;
SB = Solution.VSB;
S = [SA;SC;SB];
% Check that the sum of S components equals 1.
Ssum = simplify(collect(expand(sum(S))));
if(Ssum==1)
    fprintf('Elements of S sum up to 1\n');
end

% Accumulate initial beliefs for product A.
PPA = [PAA PA PBA];
% Accumulate initial beliefs for product B.
PPB = [PAB PB PBB];
% Compute the limiting belief for product A.
XXA = PPA * S;
% Compute the limiting belief for product A.
XXB = PPB * S;

% Simplify expressions for XA and XB given that:
% (i): PAA = PBB = 1 and
% (ii): PAB = PBA = 0.
XXA = subs(XXA,[PAA,PAB,PBA,PBB],[1 0 0 1]);
XXB = subs(XXB,[PAA,PAB,PBA,PBB],[1 0 0 1]);

% Perform some additional symbolic simplification for the temporary 
% expressions XXA and XXB.
XXA = simplify(collect(expand(XXA)));
XXB = simplify(collect(expand(XXB)));

% Important Note!!!
% Mind that at this point symbolic variables XXA and XXB hold temporary
% expressions for the limiting agents' beliefs for each product as
% functions of the network structure related parameters (TA,TB,LA,LB) and
% the initial beliefs of the consumer (PA,PB) for each product.
% These expressions will be subsequently assigned to the original symbolic
% variables for the limiting beliefs XA and XB.

% Define the demand functions for products A and B.
QA = XA - pA + M*pB - K*XB;
QB = XB - pB + M*pA - K*XA;

% Define the profit functions for the two firms A and B.
FA = QA*(pA - C) - G*TA^2;
FB = QB*(pB - C) - G*TB^2;

% Compute the first derivative of FA with respect to pA.
DFA_pA = diff(FA,pA);
% Solve the equation DFA_pA == 0 with respect to pA in order to acquire the
% Best Repsponse Function for firm A which determines the corresponding
% value pA_star as a function of pB.
% Therefore, we have that: pA = Fa(pB).
pA_star = solve(DFA_pA==0,pA);

% Compute the first derivative of FB with respect to pB.
DFB_pB = diff(FB,pB);
% Solve the equation DFB_pB == 0 with respect to pB in order to acquire the
% Best Repsponse Function for firm B which determines the corresponding
% value pB_star as a function of pA.
% Therefore, we have that pB = Fb(pA).
pB_star = solve(DFB_pB==0,pB);

% At this point we need to acquire the symbolic expressions for the optimal
% values of pA and pB (namely, pA_opt and pB_opt) by performing the following
% substitutions: 
% (i)  pA_opt may be obtained by solving with respect to pA the equation
%      pA = Fa(Fb(pA)).
% (ii) pB_opt may be obtained by solving with respect to pB the equation
%      pB = Fb(Fa(pB)).

% Make copies of the original variables in order to perform the necessary
% substitutions.
pA_star_copy = pA_star;
pB_star_copy = pB_star;

% Express variable pB in the Fa(pB) expression of pA_star as Fb(pA).
pA_star = subs(pA_star,pB,pB_star_copy);
% Solve the equation pA = Fa(Fb(pA) with respect to pA and assign the
% corresponding value to pA_opt.
pA_opt = solve(pA_star==pA,pA);

% Express variable pA in Fb(pA) as Fa(pB).
pB_star = subs(pB_star,pA,pA_star_copy);
% Solve the equation pB = Fb(Fa(pB)) with respect to pB and assign the
% corresponding value to pB_opt.
pB_opt = solve(pB_star==pB,pB);

% Substitute variables pA and pB with pA_opt and pB_opt within the 
% expressions acquired for QA and QB, thus forming the new expressions 
% QAopt and QBopt.
QA_opt = subs(QA,[pA,pB],[pA_opt,pB_opt]);
QB_opt = subs(QB,[pA,pB],[pA_opt,pB_opt]);

% Perform some further simplification on the expressions for QA_opt and
% QB_opt.
QA_opt = simplify(collect(expand(QA_opt)));
QB_opt = simplify(collect(expand(QB_opt)));

% Collect the optimal values for the prices and quantities with respect to
% XA and XB. Therefore, the quantities (pA_opt,pB_opt) and (QA_opt,QB_opt)
% will be expressed as bivariate polynomials of XA and XB.
pA_opt = collect(pA_opt,[XA,XB]);
pB_opt = collect(pB_opt,[XA,XB]);
QA_opt = collect(QA_opt,[XA,XB]);
QB_opt = collect(QB_opt,[XA,XB]);

% By performing additional symbolic operations we may write that:
% pA_opt = alpha * XA + beta * XB - gamma
% pB_opt = beta * XA + alpha * XB - gamma
% QA_opt = alpha * XA + beta * XB - gamma'
% QB_opt = beta * XA + alpha * XB - gamma'
% where the auxiliary parameters alpha, beta, gamma and gamma' are given by
%
%           K*M - 2                2*K - M                    C  
% alpha = ----------- [1] beta = ----------- [2] gamma = ----------- [3]
%           M^2 - 4                M^2 - 4                  M - 2
%
% gamma' = gamma * (M - 1)[4]
%
% Considering the differences pZ_opt - QZ_opt, where Z in {A,B}, it is easy 
% to deduce that: QZ_opt = pZ_opt - C [5] given that:
% gamma' - gamma = gamma * (M - 1) - gamma = gamma * (M - 2) = C [6].

% By setting FA_opt = FA(pA_opt,QA_opt) and FB_opt = FB(pB_opt,QB_opt) we 
% have that:
% (i):  FA_opt = (pA_opt - C)^2 - G * TA^2
% (ii): FB_opt = (pB_opt - C)^2 - G * TB^2

% Define the additional symbolic variables required.
syms alpha beta gamma


% Symbolic variables ppA and ppB which store alternative expressions for
% the quantities pA_opt and pB_opt.
ppA = alpha * XA + beta * XB - gamma;
ppB = beta * XA + alpha * XB - gamma;
% Define the symbolic expressions for the optimal values of FA_opt and
% FB_opt.
FA_opt = (ppA - C)^2 - G * TA^2;
FB_opt = (ppB - C)^2 - G * TB^2;
% Re-express symbolic variables FA_opt and FB_opt as multivariate 
% polynomials with respect to XA and XB.
FA_opt = collect(FA_opt,[XA,XB]);
FB_opt = collect(FB_opt,[XA,XB]);

% -------------------------------------------------------------------------
% SYMBOLIC COMPUTATIONS ON THE REBENUE COMPONENTS OF THE PROFIT FUNCTIONS
% -------------------------------------------------------------------------
% Obtain the optimal expressions for the revenues of the two firms.
RA_opt = (ppA - C)^2;
RB_opt = (ppB - C)^2;
% Symbolic expressions for RA_opt and RB_opt should be re-expressed as
% functions of TA and TB. The newly acquired expressions will be identified
% by Ra and Rb.
Ra = subs(RA_opt,[XA,XB],[XXA,XXB]);
Rb = subs(RB_opt,[XA,XB],[XXA,XXB]);
% Express revenue functions Ra and Rb as fractions of polynomials with
% respect to both TA and TB in the following form:
%
%        Lx
% Rx = ------ , where x in {a,b}.
%        Mx
%
[La,Ma] = numden(Ra);
[Lb,Mb] = numden(Rb);
% Check whether both fractional expressions for Ra and Rb share the same
% denumerator.
if((Ma-Mb)==0)
    fprintf('Revenue denumerators of Ra and Rb are equal\n');
end
% Extract the polynomial coefficients and corresponding monomial terms for 
% La, Ma and Lb, Mb with respect to both TA and TB.
[CLa,TLa] = coeffs(La,[TA TB]);
[CMa,TMa] = coeffs(Ma,[TA TB]);
[CLb,TLb] = coeffs(Lb,[TA TB]);
[CMb,TMb] = coeffs(Mb,[TA TB]);
% Compute the first-order derivatives of the revenue functions for the two
% firms with repsect to TA and TB accordingly.
DRa = diff(Ra,TA);
DRb = diff(Rb,TB);
% Express the above quantities as fractions of the following form:
%
%         Gx
% DRx = ------ , where x in {a,b}.
%         Sx
%
[Ga,Sa] = numden(DRa);
[Gb,Sb] = numden(DRb);

% Check whether both expressions Sa and Sb are the same.
if((Sa-Sb)==0)
    fprintf('Revenue first derivative denumerators Sa and Sb are equal\n');
end

% Extract the polynomial coefficients and corresponding monomial terms for 
% Ga, Sa and Gb, Sb with respect to both TA and TB.
[CGa,TGa] = coeffs(Ga,[TA TB]);
[CSa,TSa] = coeffs(Sa,[TA TB]);
[CGb,TGb] = coeffs(Gb,[TA TB]);
[CSb,TSb] = coeffs(Sb,[TA TB]);
% Compute the second-order derivatives of the revenue functions for the two
% firms with respect to TA and TB accoridingly.
DDRa = diff(DRa,TA);
DDRb = diff(DRb,TB);
% Express the above quantities as fractions of the following form:
%
%          Hx
% DDRx = ------ , where x in {a,b}.
%          Ox
%
[Ha,Oa] = numden(DDRa);
[Hb,Ob] = numden(DDRb);

% Check whether both expressions Sa and Sb are the same.
if((Oa-Ob)==0)
    fprintf('Revenue second derivative denumerators Oa and Ob are equal\n');
end

% Extract the polynomial coefficients and corresponding monomial terms for
% Ha, Oa and Hb, Ob with respect to both TA and TB.
[CHa,THa] = coeffs(Ha,[TA TB]);
[COa,TOa] = coeffs(Oa,[TA TB]);
[CHb,THb] = coeffs(Hb,[TA TB]);
[COb,TOb] = coeffs(Ob,[TA TB]);

% Symbolic expressions for FA_opt and FB_opt should be re-expressed as 
% functions of TA and TB. The newly acquired expressions will be identified 
% by fa and fb.
fa  = subs(FA_opt,[XA,XB],[XXA,XXB]);
fb  = subs(FB_opt,[XA,XB],[XXA,XXB]);

% -------------------------------------------------------------------------
% SYMBOLIC COMPUTATIONS FOR THE CASE LA = LB = L and PA = PB = P
% -------------------------------------------------------------------------
% Uncomment the following block of code in order to enforce LA = LB = L
% and PA = PB = P.
% Define additional symbolic variables in order to perform the addtional
% simplifications according to which: 
% LA = LB = L and
% PA = PB = P.
% syms L P
% fa = subs(fa,[LA LB PA PB],[L L P P]);
% fb = subs(fb,[LA LB PA PB],[L L P P]);

% Uncomment the following block of code in order to enforce LA = LB = 1.
% fa = subs(fa,[LA LB],[1 1]);
% fb = subs(fb,[LA LB],[1 1]);


% Extract some additional insights concerning the expressions of fa and fb.
% To this end we will express the profit functions fa and fb as fractions
% of polynomials with respect to both TA and TB as follows:
%         
%        Pa               Pb
% fa = ------  and fb = ------
%        Qa               Qb
[Pa,Qa] = numden(fa);
[Pb,Qb] = numden(fb);

% Check whether both fractional expressions for fa and fb share the same
% denumerator.
if((Qa-Qb)==0)
    fprintf('Profit denumerators Qa and Qb are equal\n');
end

% Extract the polynomial expressions and corresponding monomial terms for 
% Pa,Qa and Pb,Qb with respect to both TA and TB.
[CPa,TPa] = coeffs(Pa,[TA,TB]);
[CQa,TQa] = coeffs(Qa,[TA,TB]);
[CPb,TPb] = coeffs(Pb,[TA,TB]);
[CQb,TQb] = coeffs(Qb,[TA,TB]);

% Get the partial derivatives of fa with respect to TA and fb with respect
% to TB.
Da = diff(fa,TA);
Db = diff(fb,TB);

% Express the derivatives Da and Db as fractions of polynomial functions
% of the form:
%                   Ua             Ub  
%             Da = ----  and Db = ---- .  
%                   Va             Vb
[Ua,Va] = numden(Da);
[Ub,Vb] = numden(Db);
% Check that denumerators for both expressions (Va and Vb) are equal.
if((Va-Vb)==0)
    fprintf('Profit first derivative denumerators Va and Vb are equal!\n');
end

% Extract the polynomial expressions and corresponding monomial terms for 
% Ua,Va and Ub,Vb with respect to both TA and TB.
[CUa,TUa] = coeffs(Ua,[TA,TB]);
[CVa,TVa] = coeffs(Va,[TA,TB]);
[CUb,TUb] = coeffs(Ub,[TA,TB]);
[CVb,TVb] = coeffs(Vb,[TA,TB]);

% Given that Va = Vb reaction functions Ra(TB) and Rb(TA) for firms A and B
% may be obtained by solving:
% (i): Ua(TA,TB) = 0 which expresses TA_star as a function of TB such that
%                    TA_star = Ra(TB), and
% (ii): Ub(TA,TB) = 0 which expresses TB_star as a function of TA such that
%                     TB_star = Rb(TA).

% However, the expressions for Ua(TA,TB) and Ub(TA,TB) are proven to be 
% fourth-degree polynamials with respect to a monomial of the form: 
% TA^(d_a)*TB^(d_b) where d_a + d_b <= 4.
% Therefore, the expressions Ra and Rb are not exclusively functions of the
% variables TB and TA but they are in fact bivariate functions of TA and TB
% such that TA_star = Ra(TA,TB) and TB_star = Rb(TA,TB).
% In order to some further insight concerning the analytical forms of Ra
% and Rb we need to collect the corresponding monomial terms and associated
% coefficients for each expression.
[Ca,Ta] = coeffs(Ua,[TA,TB]);
[Cb,Tb] = coeffs(Ub,[TA,TB]);

% In order to acquire the ability to graphically represent the reaction
% curves for both utility functions Fa and Fb we need to reexpress the
% quantities Ua and Ub as polynomials of TA and TB resprectively and
% collect the corresponding coefficients and monomial terms.
[Caa,Taa] = coeffs(Ua,TA);
[Cbb,Tbb] = coeffs(Ub,TB);

% Reformulate expressions for Ua and Ub as bivariate polynomials of TA and 
% TB.
Ua = collect(Ua,[TA,TB]);
Ub = collect(Ub,[TA,TB]);

% Since Ua(TA,TB) = 0 and Ub(TA,TB) = 0, we may use an additional equation
% such that Uab = Ua(TA,TB) - Ub(TA,TB) = 0.
Uab = Ua - Ub;
% Get the coefficients and corresponding monomial terms of the polynomial
% Uab(TA,TB).
[Cab,Tab] = coeffs(Uab,[TA,TB]);
% Simplify the expression Uab.
Uab = simplify(collect(expand(Uab)));
% Factorize the expression Uab with respect to TA and TB.
Fab = factor(Uab,[TA,TB]);

% Compute the second derivative of fa and fb with respect to TA and TB
% respectivelly.
DDa = diff(Da,TA);
DDb = diff(Db,TB);

% Rexpress quantities DDa and DDb as fractions of the following form:
%         Wa                 Wb
% DDa = ------   and DDb = ------
%         Za                 Zb
[Wa,Za] = numden(DDa);
[Wb,Zb] = numden(DDb);

% Check that denumerators for both expressions (Za and Zb) are equal.
if((Za-Zb)==0)
    fprintf('Profit second derivative denumerators Za and Zb are equal!\n');
end

% Extract the polynomial expressions and corresponding monomial terms for 
% Wa,Za and Wb,Zb with respect to both TA and TB.
[CWa,TWa] = coeffs(Wa,[TA,TB]);
[CZa,TZa] = coeffs(Za,[TA,TB]);
[CWb,TWb] = coeffs(Wb,[TA,TB]);
[CZb,TZb] = coeffs(Zb,[TA,TB]);