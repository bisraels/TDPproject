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
