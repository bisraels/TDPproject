% Try calculating correlation functions with matricies

K = [-25,  21,     31,         0,     0,       0,       0;...   %1
    12,    -93,    2852/91,    42,    52,      0,       0;...   %2
    13,    23,  -(31+2852/91), 0,     0,       0,       0;...   %3
    0,     24,  0,             -42,   0,       0,       0;...   %4
    0,     25,  0,              0,    -108,    65,      0;...   %5
    0,     0,   0,              0,      56,  -(65+67),  76;...  %6
    0,     0,   0,              0,       0,     67,     -(76);] %7
A = [.1:.1:.7];
N = length(K);


[Vec, lambda] = eig(K);

elapsedTime = toc;
% if clockMode == 1
%     task_str = ['solve the ' num2str(N) 'N-state eigensystem.'];
%     disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
% end

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);  % This is our unitary transform matrix.

V = VecSorted;

syms t
expLam = diag(vpa(exp(diag(lambdaSorted) * t)));

initCond = eye(7);

p = vpa(V * expLam * inv(V) * initCond);

peq = diag(p);

A1 = 0.8;
A2 = 0.5;
A3 = 0.35;
A4 = 0.3;
A5 = 0.4;
A6 = 0.45;
A7 = 0.2;
A = diag([A1; A2; A3; A4; A5; A6; A7]);

c2 = sum(A * p * A * peq);


%%
syms tau1 tau2 tau3
t = tau1;
p_tau1 = subs(p);  % Same as for tau 3

t = tau2;
p_tau2 = subs(p);

tau3 = tau1;
p_tau3 = subs(p);


t = 300;
Peq = double(subs(peq));
% 
% 
c4 = A * p_tau3 * A * p_tau2 * A * p_tau1 * A * Peq;
c4_sum = sum(A * p_tau3 * A * p_tau2 * A * p_tau1 * A * Peq);
% 
% tau2 = 0;
% tau1 = 1:100;
% tau3 = tau1';
% 
% C4 = double(subs(c4));
% size(C4)
% 
% 
% C4_mat = C4' * C4;
% surf(C4_mat);
% set




%% HAVE NOT FIGURED OUT HOW TO MAKE THIS WORK YET

tau1 = 1:100;
tau2 = 1;
tau3 = tau1;
P_tau1 = k2P(K,tau1);
P_tau2 = k2P(K,tau2);
P_tau3 = k2P(K,tau3);

Peq = diag(Peq);

P_tau3_p = permute(P_tau3,[3,1,2]) ;

C2_matSum = {};
C2_matSum_temp = [];
for j = 1:length(tau1)
for i = 1:length(tau1)
    C2_matSum_temp = reshape(P_tau3_p(i,:,:),N,N) * A * A * A * A * P_tau2 * Peq * reshape(P_tau1(:,:,j),N,N);
    C2_matSum(i,j) = C2_matSum_temp;
end
end
size(C2_matSum)

% P_tau3_rs = reshape(P_tau3, 100, 7, 7);
% 
% P_tau3_rs2 = reshape(P_tau3(:),100,7,7);

%%

P_prod = P_tau3_p * P_tau1











