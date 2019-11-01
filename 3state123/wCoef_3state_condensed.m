clockMode = 1;

% Define matrix of c's for our constants
c = sym('c', [3,3]);

% Define matrix of the eigenvector components
V = sym('v',[3,3]).' ;

% Use identity matrix to define the possible conditions
cond = eye(3);

% Use these to construct our equations
eq1 = cond(:,1) == V * c(:,1) ;
eq2 = cond(:,2) == V * c(:,2) ;
eq3 = cond(:,3) == V * c(:,3) ;

% Define the eqilibrium probabilities
P1_eq = (c(1,1) * V(1,1)); % from P11
P2_eq = (c(1,2) * V(2,1)); % from P22
P3_eq = (c(1,3) * V(3,1)); % from P33


% Additional constraint such that Peq sums to 1.
eq_sum = P1_eq + P2_eq + P3_eq == 1;

% eqs = [eq1(1),eq1(2),eq1(3),...
%     eq2(1),eq2(2),eq2(3),...
%     eq3(1),eq3(2),eq3(3),...
%     eq_sum];
eqs = [eq1, eq2, eq3, eq_sum];

% Build list of vars - writing out by hand is prone to  mistakes
vars = c(:);

% if length(vars) == (length(c)*length(c))
%     disp('We have all the variables!')
% else
%     disp('We have missing/extra variables.')
% end
    
% Solve eqns for vars
tic
sols = solve(eqs, vars)
toc


% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * V(1,1)); 
P2_eq = (sols.c1_2 * V(2,1)); 
P3_eq = (sols.c1_3 * V(3,1)); 

c1_1 = sols.c1_1;
c2_1 = sols.c2_1;
c3_1 = sols.c3_1;

c1_2 = sols.c1_2;
c2_2 = sols.c2_2;
c3_2 = sols.c3_2;

c1_3 = sols.c1_3;
c2_3 = sols.c2_3;
c3_3 = sols.c3_3;

save('wCoef_3state_condensed.mat','sols',...
    'c1_1', 'c1_2', 'c1_3',...
    'c2_1', 'c2_2', 'c2_3',...
    'c3_1', 'c3_2', 'c3_3',...
    'P1_eq', 'P2_eq','P3_eq');

%%
load('wCoef_3state_condensed.mat',...
    'c1_1', 'c1_2', 'c1_3',...
    'c2_1', 'c2_2', 'c2_3',...
    'c3_1', 'c3_2', 'c3_3',...
    'P1_eq', 'P2_eq','P3_eq');
% Test the solutions with calcuation of 2 point TCF
%--------------------------------------------------------------------------
% Pick random values for the rate constants and FRET states
%--------------------------------------------------------------------------
[k12,k13,k21,k23,k31,A1,A2,A3,k32,time] = paramSim_3state123_cyclical();

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors
%--------------------------------------------------------------------------
K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];


% if clockMode == 1
%     tic
% end
[Vec, lambda] = eig(K);   
% %Display the amount of time a process took. Begins at the last tic.
% if clockMode == 1
%     elapsedTime = toc;
%     task_str = 'Calculate the eigenvalues';
%     disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
% end

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

lam2 = double(lambdaSorted(2,2));
lam3 = double(lambdaSorted(3,3));

evec1 = VecSorted(:,1);%We dont need this eigenvector
evec2 = VecSorted(:,2);
evec3 = VecSorted(:,3);


v1_1 = evec1(1);
v1_2 = evec1(2);
v1_3 = evec1(3);
v2_1 = evec2(1);
v2_2 = evec2(2);
v2_3 = evec2(3);
v3_1 = evec3(1);
v3_2 = evec3(2);
v3_3 = evec3(3);

% Calculate the expansion coefficients in terms of the eigevcetor
% components
c1_1 = subs(c1_1);

c2_1 = subs(c2_1);
c3_1 = subs(c3_1);

c1_2 = subs(c1_2);
c2_2 = subs(c2_2);
c3_2 = subs(c3_2);

c1_3 = subs(c1_3);
c2_3 = subs(c2_3);
c3_3 = subs(c3_3);

P1_eq = subs(P1_eq);
P2_eq = subs(P2_eq);
P3_eq = subs(P3_eq);

% Calculate the NxN conditional Probabilities
p1_1 = P1_eq + c2_1*v2_1*exp(-lam1*time) + c3_1*v3_1*exp(-lam2*time);
p1_2 = P1_eq + c2_2*v2_1*exp(-lam1*time) + c3_2*v3_1*exp(-lam2*time);
p1_3 = P1_eq + c2_3*v2_1*exp(-lam1*time) + c3_3*v3_1*exp(-lam2*time);
p2_1 = P2_eq + c2_1*v2_2*exp(-lam1*time) + c3_1*v3_2*exp(-lam2*time);
p2_2 = P2_eq + c2_2*v2_2*exp(-lam1*time) + c3_2*v3_2*exp(-lam2*time);
p2_3 = P2_eq + c2_3*v2_2*exp(-lam1*time) + c3_3*v3_2*exp(-lam2*time);
p3_1 = P3_eq + c2_1*v2_3*exp(-lam1*time) + c3_1*v3_3*exp(-lam2*time);
p3_2 = P3_eq + c2_2*v2_3*exp(-lam1*time) + c3_2*v3_3*exp(-lam2*time);
p3_3 = P3_eq + c2_3*v2_3*exp(-lam1*time) + c3_3*v3_3*exp(-lam2*time);


% Compute the 2 point TCF (Using conditinal Probabilities)
% Subtract mean values
Amean = P1_eq*A1 + P2_eq*A2 + P3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

% Calculate C2_sim
tic
C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
toc
