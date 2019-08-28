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
k12 = 44;
k13 = 0;
k20 = 5;
k21 = 60;
k23 = 7;
k30 = 0;
k31 = 0;
k32 = 87;
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

P0EQ = vpa(subs(P0eq));

P10(t) = vpa(subs(P10));
P11(t) = vpa(subs(P11));
P12(t) = vpa(subs(P12));
P13(t) = vpa(subs(P13));

P1EQ = vpa(subs(P1eq));

P20(t) = vpa(subs(P20));
P21(t) = vpa(subs(P21));
P22(t) = vpa(subs(P22));
P23(t) = vpa(subs(P23));

P2EQ = vpa(subs(P2eq));

P30(t) = vpa(subs(P30));
P31(t) = vpa(subs(P31));
P32(t) = vpa(subs(P32));
P33(t) = vpa(subs(P33));

P3EQ = vpa(subs(P3eq));

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

%--------------------------------------------------------------------------
% Two point TCF:
%--------------------------------------------------------------------------
% syms A0 A1 A2 A3
syms t

% This is a vector for the FRET values - assign values
A0 = 0.34;  % These are just guesses for now
A1 = 0.42;
A2 = 0.51;
A3 = 0.67;

A = [A0; A1; A2; A3];

% Matrix of conditional probabilities Pi-->j with i is initial condition
% cP = [P00(t), P01(t), P02(t), P03(t);...
%      P10(t), P11(t), P12(t), P13(t);...
%      P20(t), P21(t), P22(t), P23(t);...
%      P30(t), P31(t), P32(t), P33(t)];
 cP = [P00(t), P10(t), P20(t), P30(t);...
       P01(t), P11(t), P21(t), P31(t);...
       P02(t), P12(t), P22(t), P32(t);...
       P03(t), P13(t), P23(t), P33(t)];
% Row = final state & Column = initial condition

%%
load('sym_prob.mat')
% This gives the equilibrium probabilties
% for this state we are using:
% Pi_Cn4 where i = 0, 1, 2, 3
 
% P0_eq = vpa(subs(P0_Cn4));
% P1_eq = vpa(subs(P1_Cn4));
% P2_eq = vpa(subs(P2_Cn4));
% P3_eq = vpa(subs(P3_Cn4));
% 
% Peq = [P0_eq; P1_eq; P2_eq; P3_eq];

P0EQ = P00(inf);
P1EQ = P11(inf);
P2EQ = P22(inf);
P2EQ = P33(inf);

Peq = [P0EQ; P1EQ; P2EQ; P3EQ];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

%%
Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


C2(t) = 0;
for i = 1:numel(A)
    for j = 1:numel(A)
        
        C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
        C2(t) = C2(t) + C2temp(t);
    end
end
%%

if double(C2(0)) == double(msq)
    disp('Mean of the square <A^2> matches C2(t=0)!')
else 
    disp('Problem: mean of the square <A^2> DOES NOT match C2(t=0)')
end

if C2(10^20) == sqm
    disp('Square of the mean <A>^2 matches C2(t=BIG)!')
else 
    disp('Problem: square of the mean <A>^2 DOES NOT match C2(t=0)')
    
end

%%
%--------------------------------------------------------------------------
% Four point TCF:
%--------------------------------------------------------------------------


% All same info from above

% C4(t) = 0;
% 
% for l = 1:numel(A)
%     for k = 1:numel(A)
%         for j = 1:numel(A)
%             for i = 1:numel(A)
%                 
%                 C4temp(t) = A(l) * cP(l,k) * A(k) * cP(k,j) * A(j) * cP(j,i) * A(i) * Peq(i);
%                 C4(t) = C4(t) + C4temp(t);
%             end
%         end
%     end
% end

%%
% Plot TCF's
close all

figure
% subplot(1,2,1)
fplot(C2(t),[0,1])
title('Two point TCF')

% hold on
% subplot(1,2,2)
% fplot(C4(t),[0,1])
% title('Four point TCF')
% hold off


