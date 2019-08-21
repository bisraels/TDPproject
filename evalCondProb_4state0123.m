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



