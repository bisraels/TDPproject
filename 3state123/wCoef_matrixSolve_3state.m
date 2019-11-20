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

saveMode = 0;

N = 3;

% Define matrix of c's for our constants
c = sym('c', [N,N]);

% Define matrix of the eigenvector components
v = sym('v',[N,N]).' ;

% Use identity matrix to define the possible conditions
P0 = eye(N);

% Find the inverse of the matrix of eigenvectors
tic
v_inv = inv(v);
elapsedTime = toc;
disp(['Time to find inverse of V = ' num2str(elapsedTime)]);

tic
% Find the matrix of all possible initial conditions
% c_mat  = v_inv * P0;    %Slower
c_mat  = P0 \ v;          %Faster
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
if saveMode ==1 
save('wcoef_matSolve_3state_cMatrix.mat','c_mat')
end
%--------------------------------------------------------------------------
% Pick random values for the rate constants and FRET states
%--------------------------------------------------------------------------
% [k12,k13,k21,k23,k31,A1,A2,A3,k32,~] = paramSim_3state123_cyclical();

k12 = 12; k13 = 13; k21 = 21; k31 = 31; k23 = 23;
k32 = k12*k23*k31/(k13*k21);

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors of K matrix
%--------------------------------------------------------------------------

K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];


[evec, lam] = eig(K);   % Need to sort these to test this solution properly.

[lam_sort,ind] = sort(diag(lam),'descend');
evec_sort = evec(:,ind);

evec_1 = evec_sort(:,1);
evec_2 = evec_sort(:,2);
evec_3 = evec_sort(:,3);


v1_1 = evec_1(1);
v1_2 = evec_1(2);
v1_3 = evec_1(3);

v2_1 = evec_2(1);
v2_2 = evec_2(2);
v2_3 = evec_2(3);

v3_1 = evec_3(1);
v3_2 = evec_3(2);
v3_3 = evec_3(3);

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
cPj_1(t) = vpa(subs(v * exp_LamT * v_inv * P0(:,1)));
% Prob: 2  --> j
cPj_2(t) = vpa(subs(v * exp_LamT * v_inv * P0(:,2)));
% Prob: 3  --> j
cPj_3(t) = vpa(subs(v * exp_LamT * v_inv * P0(:,3)));

elapsedTime = toc;
disp(['Time to evaluate cPj_i function of t = ' num2str(elapsedTime)]);



P = vpa(subs(v * exp_LamT * v_inv * P0));