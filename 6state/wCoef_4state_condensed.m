
N = 4;

% Define matrix of c's for our constants
c = sym('c', [N,N]);

% Define matrix of the eigenvector components
V = sym('v',[N,N]).' ;

% Use identity matrix to define the possible conditions
cond = eye(N);

% Use these to construct our equations
eq1 = cond(:,1) == V * c(:,1) ;
eq2 = cond(:,2) == V * c(:,2) ;
eq3 = cond(:,3) == V * c(:,3) ;
eq4 = cond(:,4) == V * c(:,4) ;


% Define the eqilibrium probabilities
P1_eq = (c(1,1) * V(1,1)); % from P11
P2_eq = (c(1,2) * V(2,1)); % from P22
P3_eq = (c(1,3) * V(3,1)); % from P33
P4_eq = (c(1,4) * V(4,1)); % from P44

% Additional constraint such that Peq sums to 1.
eq_sum = P1_eq + P2_eq + P3_eq + P4_eq == 1;

eqs = [eq1, eq2, eq3, eq4, eq_sum];

% Build list of vars - writing out by hand is prone to  mistakes
vars = c(:);


if length(vars) == (length(c)*length(c))
    disp('We have all the variables!')
else
    disp('We have missing/extra variables.')
end
    
    
% Solve eqns for vars
tic
sols = solve(eqs, vars)
toc
%Takes about .18 seconds

% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * V(1,1)); 
P2_eq = (sols.c1_2 * V(2,1)); 
P3_eq = (sols.c1_3 * V(3,1)); 
P4_eq = (sols.c1_4 * V(4,1)); 

c1_1 = sols.c1_1;
c2_1 = sols.c2_1;
c3_1 = sols.c3_1;
c4_1 = sols.c4_1;

c1_2 = sols.c1_2;
c2_2 = sols.c2_2;
c3_2 = sols.c3_2;
c4_2 = sols.c4_2;

c1_3 = sols.c1_3;
c2_3 = sols.c2_3;
c3_3 = sols.c3_3;
c4_3 = sols.c4_3;

c1_4 = sols.c1_4;
c2_4 = sols.c2_4;
c3_4 = sols.c3_4;
c4_4 = sols.c4_4;


save('wCoef_4state_condensed.mat','sols',...
    'c1_1', 'c1_2', 'c1_3', 'c1_4', ...
    'c2_1', 'c2_2', 'c2_3', 'c2_4', ...
    'c3_1', 'c3_2', 'c3_3', 'c3_4', ...
    'c4_1', 'c4_2', 'c4_3', 'c4_4', ...
    'P1_eq', 'P2_eq','P3_eq','P4_eq');
