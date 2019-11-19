%__________________________________________________________________________
% AUTHOR:   Brett Israels and Claire Albrecht
%
% CREATED:  November 2019
%
% PURPOSE:  Solve the master equation for a rate matrix K numerically
%
% INPUT     k (rate matrix)
%
% OUTPUT:   P (matrix of conditional probabilites)
%           P(j,i,t) is the probability of going from i -> j in time t
%
% MODLOG:  BI 20191112 Inverse done with '\' not inv() (faster)
%          CA
%__________________________________________________________________________

function [P] = k2P(K,time)
switch nargin
    case 0
        %                 Simulate a K matrix
        %                 3 state
        %                 k12 = 12; k13 = 13; k21 = 21; k31 = 31; k23 = 23;
        %                 k32 = k12*k23*k31/(k13*k21);
        %
        %                 K = [(-k12 - k13), k21, k31;...
        %                     k12, (-k21 - k23 ), k32;...
        %                     k13, k23, (-k31-k32);];
        %
        %                 8state
        K = [...
            -(12 + 13), 21, 31, 0, 0, 0, 0, 0;...
            12, -(21 + 23 + 24 + 25), (12*23*31/(13*21)), 42, 52, 0, 0, 0;...
            13, 23, -(31 + ((12*23*31/(13*21)))), 0, 0, 0, 0, 0;...
            0, 24, 0, -42, 0, 0, 0, 0;...
            0, 25, 0, 0, -(52 + 56), 65, 0, 0;...
            0, 0, 0, 0, 56, -(65 + 67 + 68), 76, 86;...
            0, 0, 0, 0, 0, 67, -(76+78), 87;...
            0, 0, 0, 0, 0, 68, 78, -(86 + 87);...
            ];
        
        Npts = 150;
        time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;
        
    case 1 %If you give it 1 arguement , assume the gp32 = 0.
        time = 0;
end
% Determine the number of states in the system
N = length(K);

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors of K matrix
%--------------------------------------------------------------------------
[Evec, Lam_unsorted] = eig(K);

[Lam,ind] = sort(diag(Lam_unsorted),'descend'); %Lam is numerical eigenvalue
V = Evec(:,ind);        % V is the mode matrix (Each column is an eigenvector)

% Find the inverse of the matrix of eigenvectors
tic
V_inv = inv(V);
elapsedTime = toc;
disp(['Time to find inverse of V = ' num2str(elapsedTime)]);

% Use identity matrix to define the Initial conditions (Boundary conditions)
InitCond = eye(N);

% Find the matrix of all possible initial conditions
tic
C  = V \ InitCond; %V_inv * InitCond;
elapsedTime = toc;
disp(['Time to find C matrix = ' num2str(elapsedTime)]);
% Format:
%      c_mat = c1_1, c1_2, ..., c1_n;
%              c2_1, c2_2, ..., c2_n;
%               .    .           .
%               .      .         .
%               .        .       .
%              cn_1     ...     cn_n
%

%--------------------------------------------------------------------------
% Solve for the conditional probabilities
%--------------------------------------------------------------------------
syms t

exp_LamT = diag(exp(Lam * t));
% tic
% cP(t) = vpa(subs(V * exp_LamT * V_inv * InitCond));
% elapsedTime = toc;
% disp(['Time to evaluate cP matrix as function of t = ' num2str(elapsedTime)]);


tic
%equivalent to P = U * exp(lam*t) * C
p = vpa(subs(V * exp_LamT * V_inv * InitCond));

%P(j,i) is the probability of going from i --> j in time t
elapsedTime = toc;
disp(['Time to evaluate P(j,i) function of t = ' num2str(elapsedTime)]);

% %Note: To substitute in time, do so like this:
t = time;
% t = 0:.1:1;
% tic
P = reshape(subs(p(:)),[N N numel(t)]);
% toc