% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE: Evaluate the four state conditional probabilities with set of
% rates
%
% MODIFICATIONS:
%--------------------------------------------------------------------------

%load('symCondProb_4state0123.mat')        % This is from the old code
load('symCondProb_4state0123_fixed.mat')   % This is the output of the fixed code.

syms t

k01 = 1;
%k02 = 2;
k03 = 0;
k10 = 3;
k12 = 4;
k13 = 0;
k20 = 5;
k21 = 6;
k23 = 7;
k30 = 0;
k31 = 0;
k32 = 8;
% NOTE: I put in the rates we know are zero in by hand.

% Detailed balance condition
% k20 = (k10 * k02 * k21)/(k10 * k12);
k02 = (k01 * k12 * k20)/(k10 * k21);

% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.
P00(t) = vpa(subs(P00));    
P01(t) = vpa(subs(P01));
P02(t) = vpa(subs(P02));
P03(t) = vpa(subs(P03));

P10(t) = vpa(subs(P10));
P11(t) = vpa(subs(P11));
P12(t) = vpa(subs(P12));
P13(t) = vpa(subs(P13));

P20(t) = vpa(subs(P20));
P21(t) = vpa(subs(P21));
P22(t) = vpa(subs(P22));
P23(t) = vpa(subs(P23));

P30(t) = vpa(subs(P30));
P31(t) = vpa(subs(P31));
P32(t) = vpa(subs(P32));
P33(t) = vpa(subs(P33));

% NOTE: When using these expressions, when you plug in a value of t, the
% expression will still be a of class 'sym'.
% To evaluate these expressions as doubles use the following:
%       If t = 0:
%       P00_t1 = double(P00(0))
% Now P00_t1 will be a double


%%
% Calculate TCF 
% Need:
%       - Pji(t) calcualted above
%       - Pi_eq equilibruim populations (from TDP_probcalc)
%       - Ai values of FRET states

% Two point TCF:
% syms A0 A1 A2 A3

% This is a vector for the FRET values - assign values
A0 = 0.34;
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;

A = [A0, A1, A2, A3];
% Matrix of conditional probabilities
cP = [P00(t), P01(t), P02(t), P03(t);...
     P10(t), P11(t), P12(t), P13(t);...
     P20(t), P21(t), P22(t), P23(t);...
     P30(t), P31(t), P32(t), P33(t)];

load('sym_prob.mat')
% This gives the equilibrium probabilties
% for this state we are using:
% Pi_Cn4 where i = 0, 1, 2, 3
 
P0_eq = vpa(subs(P0_Cn4));
P1_eq = vpa(subs(P1_Cn4));
P2_eq = vpa(subs(P2_Cn4));
P3_eq = vpa(subs(P3_Cn4));

Peq = [P0_eq; P1_eq; P2_eq; P3_eq];

C2(t) = 0;
for i = 1:numel(A)
    for j = 1:numel(A)
        
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2(t) = C2(t) + C2temp(t);
    end
end


