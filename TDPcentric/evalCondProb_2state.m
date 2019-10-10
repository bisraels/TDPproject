% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: August 2019
%
% PURPOSE: Evaluate the two state conditional probabilities with set of
% rates
%
% MODIFICATIONS:


load('symCondProb_2state')

syms t

k12 = 1/4;
k21 = 1/2;

P11(t) = subs(P11);
P12(t) = subs(P12);
P21(t) = subs(P21);
P22(t) = subs(P22);

% Calculate two point TCF for the two state system

load('sym_prob.mat')

P1_eq = vpa(subs(P1_n2));
P2_eq = vpa(subs(P2_n2));

Peq = [P1_eq; P2_eq];

if sum(Peq) == 1
    disp('Equilibrium probabilities sum to 1!')
else
    disp('Problem: Equilibrium probabilities DO NOT sum to 1.')
end

A1 = 0.32;
A2 = 0.47;
A = [A1; A2];

Amean = sum(A.*Peq);
A = A - Amean;

msq = sum((A.^2).*Peq);     % square of mean <A^2>
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


cP = [P11(t), P21(t);...
      P12(t), P22(t)];

  C2(t) = 0;
  for i = 1:numel(A)
      for j = 1:numel(A)
          C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
          C2(t) = C2(t) + C2temp(t);
      end
  end
  
  
if double(C2(0)) == double(msq)
    disp('Mean of the square < A^2 > matches C2(t=0)!')
else 
    disp('Problem: mean of the square < A^2 > DOES NOT match C2(t=0)')
end

if C2(10^20) == sqm
    disp('Square of the mean < A >^2 matches C2(t=BIG)!')
else 
    disp('Problem: square of the mean < A >^2 DOES NOT match C2(t=BIG)')
end

fplot(C2(t),[0 10])
