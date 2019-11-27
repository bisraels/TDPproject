%  function simCnP(N)

N=3;

% Construct symbolic matrix of k's
k = sym('k%d%d',[N,N]).';
for i = 1:N
    k(i,i) = 0;
    k(i,i) = -sum(k(:,i));
end


% Define parameters 
[k12,k13,k21,k23,k31,A1,A2,A3,k32,time] = paramSim_3state123_cyclical();


% Substitute values in for K
tic
K = subs(k);
toc

% Calculate conditional probability matrix and modal matrix
[P, V, time] = k2P(K,time);

% Calculate C2 and C4
[time, C2, C4] = P2C(P, K, time);













% Rates  for  8 state
% k12 = 12;
% k13 = 13;
% k21 = 21;
% k23 = 23;
% k24 = 24;
% k25 = 25;
% k31 = 31;
% k42 = 42;
% k52 = 52;
% k56 = 56;
% k65 = 65;
% k67 = 67;
% k68 = 68;
% k76 = 76;
% k78 = 78;
% k86 = 86;
% k87 = 87;
