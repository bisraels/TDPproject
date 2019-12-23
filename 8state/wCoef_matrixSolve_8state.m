%__________________________________________________________________________
% AUTHOR:   Claire Albrecht
%
% CREATED:  November 2019
%
% PURPOSE:  Solve for weighting coefficients or conditional probabilities
%           using matrices.
%
% INPUT     None
%
% OUTPUT:   
%
% MODIFICATION LOG:
%__________________________________________________________________________


N = 8;

% Define matrix of c's for our constants
c = sym('c', [N,N]);

% Define matrix of the eigenvector components
V = sym('v',[N,N]).' ;

% Use identity matrix to define the possible conditions
cond = eye(N);

% Find the inverse of the matrix of eigenvectors
tic
V_inv = inv(V);
elapsedTime = toc;
disp(['Time to find inverse of V = ' num2str(elapsedTime)]);


% Find the matrix of all possible initial conditions
c_mat  = V_inv * cond;
elapsedTime = toc;
disp(['Time to find c matrix = ' num2str(elapsedTime)]);
% Format: 
%      c_mat = c1_1, c1_2, ..., c1_n;
%              c2_1, c2_2, ..., c2_n;
%               .    .           .
%               .      .         .
%               .        .       .
%              cn_1     ...     cn_n
%
% Select by COLUMNS for each equation.
save('wcoef_matSolve_8state_cMatrix.mat','c_mat')

return

%% Check results with arbitrary matrix

K = [-(12 + 13), 21, 31, 0, 0, 0, 0;...
    12, -(21 + 23 + 24 + 25), (12*23*31/(13*21)), 42, 52, 0, 0;...
    13, 23, -(31 + ((12*23*31/(13*21)))), 0, 0, 0, 0;...
    0, 24, 0, -42, 0, 0, 0;... 
    0, 25, 0, 0, -(52 + 56), 65, 0;...
    0, 0, 0, 0, 56, -(65 + 67), 76;...
    0, 0, 0, 0, 0, 67, -76];

[evec, lam] = eig(K);   % Need to sort these to test this solution properly.


[lam_sort,ind] = sort(diag(lam),'descend');
evec_sort = evec(:,ind);

evec_1 = evec_sort(:,1);
evec_2 = evec_sort(:,2);
evec_3 = evec_sort(:,3);
evec_4 = evec_sort(:,4);
evec_5 = evec_sort(:,5);
evec_6 = evec_sort(:,6);
evec_7 = evec_sort(:,7);


v1_1 = evec_1(1);
v1_2 = evec_1(2);
v1_3 = evec_1(3);
v1_4 = evec_1(4);
v1_5 = evec_1(5);
v1_6 = evec_1(6);
v1_7 = evec_1(7);

v2_1 = evec_2(1);
v2_2 = evec_2(2);
v2_3 = evec_2(3);
v2_4 = evec_2(4);
v2_5 = evec_2(5);
v2_6 = evec_2(6);
v2_7 = evec_2(7);

v3_1 = evec_3(1);
v3_2 = evec_3(2);
v3_3 = evec_3(3);
v3_4 = evec_3(4);
v3_5 = evec_3(5);
v3_6 = evec_3(6);
v3_7 = evec_3(7);

v4_1 = evec_4(1);
v4_2 = evec_4(2);
v4_3 = evec_4(3);
v4_4 = evec_4(4);
v4_5 = evec_4(5);
v4_6 = evec_4(6);
v4_7 = evec_4(7);

v5_1 = evec_5(1);
v5_2 = evec_5(2);
v5_3 = evec_5(3);
v5_4 = evec_5(4);
v5_5 = evec_5(5);
v5_6 = evec_5(6);
v5_7 = evec_5(7);

v6_1 = evec_6(1);
v6_2 = evec_6(2);
v6_3 = evec_6(3);
v6_4 = evec_6(4);
v6_5 = evec_6(5);
v6_6 = evec_6(6);
v6_7 = evec_6(7);

v7_1 = evec_7(1);
v7_2 = evec_7(2);
v7_3 = evec_7(3);
v7_4 = evec_7(4);
v7_5 = evec_7(5);
v7_6 = evec_7(6);
v7_7 = evec_7(7);


% Evaluate our matrix of C's in terms of the defined values.
tic 
C_mat = double(subs(c_mat));
elapsedTime = toc;
disp(['Time to evaluate C matrix = ' num2str(elapsedTime)]);




%%
% Can also solve for cP(t) in this method.

syms t

exp_LamT = diag(exp(lam_sort * t));
% tic 
% cP(t) = vpa(subs(V * exp_LamT * V_inv * cond));
% elapsedTime = toc;
% disp(['Time to evaluate cP matrix as function of t = ' num2str(elapsedTime)]);

tic

% Prob: 1  --> j
cPj_1(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,1)));
% Prob: 2  --> j
cPj_2(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,2)));
% Prob: 3  --> j
cPj_3(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,3)));
% Prob: 4  --> j
cPj_4(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,4)));
% Prob: 5  --> j
cPj_5(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,5)));
% Prob: 6  --> j
cPj_6(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,6)));
% Prob: 7  --> j
cPj_7(t) = vpa(subs(V * exp_LamT * V_inv * cond(:,7)));

elapsedTime = toc;
disp(['Time to evaluate cPj_i function of t = ' num2str(elapsedTime)]);


