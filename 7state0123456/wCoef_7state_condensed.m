% Define matrix of c's for our constants
c = sym('c', [7,7]);

% Define matrix of the eigenvector components
V = sym('v',[7,7]).' ;

% Use identity matrix to define the possible conditions
cond = eye(7);

% Use these to construct our equations
eq1 = cond(:,1) == V * c(:,1) ;
eq2 = cond(:,2) == V * c(:,2) ;
eq3 = cond(:,3) == V * c(:,3) ;
eq4 = cond(:,4) == V * c(:,4) ;
eq5 = cond(:,5) == V * c(:,5) ;
eq6 = cond(:,6) == V * c(:,6) ;
eq7 = cond(:,7) == V * c(:,7) ;

% Define the eqilibrium probabilities
% P1_eq = (c1_1 * v1_1); % from P11
% P2_eq = (c1_2 * v1_2); % from P22
% P3_eq = (c1_3 * v1_3); % from P33
% P4_eq = (c1_4 * v1_4); % from P44
% P5_eq = (c1_5 * v1_5); % from P55
% P6_eq = (c1_6 * v1_6); % from P66
% P7_eq = (c1_7 * v1_7); % from P77

P1_eq = (c(1,1) * V(1,1)); % from P11
P2_eq = (c(1,2) * V(2,1)); % from P22
P3_eq = (c(1,3) * V(3,1)); % from P33
P4_eq = (c(1,4) * V(4,1)); % from P44
P5_eq = (c(1,5) * V(5,1)); % from P55
P6_eq = (c(1,6) * V(6,1)); % from P66
P7_eq = (c(1,7) * V(7,1)); % from P77


% Additional constraint such that Peq sums to 1.
eq_sum = P1_eq + P2_eq + P3_eq + P4_eq + P5_eq + P6_eq + P7_eq == 1;

eqs = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq_sum];
% vars  = [c1_1, c2_1, c3_1, c4_1, c5_1, c6_1, c7_1, c1_2, c2_2, c3_2, c4_2, c5_2, c6_2, c7_2, c1_3, c2_3, c3_3, c4_3, c5_3, c6_3, c7_3, c1_4, c2_4, c3_4, c4_4, c5_4, c6_4, c7_4, c1_5, c2_5, c3_5, c4_5, c5_5, c6_5, c7_5, c1_6, c2_6, c3_6, c4_6, c5_6, c6_6, c7_6, v1_7, v2_7, v3_7, v4_7, v5_7, v6_7, v7_7];


% Build list of vars - writing out by hand is prone to  mistakes
vars = sym.empty();
for j = 1:length(c)
    for i = 1:length(c)
        temp(i) = c(i,j);
        vars_temp(i) = temp(i);
    end
    vars = [vars, vars_temp];
end

if length(vars) == (length(c)*length(c))
    disp('We have all the variables!')
else
    disp('We have missing/extra variables.')
end
    
    
% Solve eqns for vars
tic
sols = solve(eqs, vars)
toc

% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * V(1,1)); 
P2_eq = (sols.c1_2 * V(2,1)); 
P3_eq = (sols.c1_3 * V(3,1)); 
P4_eq = (sols.c1_4 * V(4,1)); 
p5_eq = (sols.c1_5 * V(5,1));
P6_eq = (sols.c1_6 * V(6,1)); 
P7_eq = (sols.c1_7 * V(7,1)); 

c1_1 = sols.c1_1;
c2_1 = sols.c2_1;
c3_1 = sols.c3_1;
c4_1 = sols.c4_1;
c5_1 = sols.c5_1;
c6_1 = sols.c6_1;
c7_1 = sols.c7_1;

c1_2 = sols.c1_2;
c2_2 = sols.c2_2;
c3_2 = sols.c3_2;
c4_2 = sols.c4_2;
c5_2 = sols.c5_2;
c6_2 = sols.c6_2;
c7_2 = sols.c7_2;

c1_3 = sols.c1_3;
c2_3 = sols.c2_3;
c3_3 = sols.c3_3;
c4_3 = sols.c4_3;
c5_3 = sols.c5_3;
c6_3 = sols.c6_3;
c7_3 = sols.c7_3;

c1_4 = sols.c1_4;
c2_4 = sols.c2_4;
c3_4 = sols.c3_4;
c4_4 = sols.c4_4;
c5_4 = sols.c5_4;
c6_4 = sols.c6_4;
c7_4 = sols.c7_4;

c1_5 = sols.c1_5;
c2_5 = sols.c2_5;
c3_5 = sols.c3_5;
c4_5 = sols.c4_5;
c5_5 = sols.c5_5;
c6_5 = sols.c6_5;
c7_5 = sols.c7_5;

c1_6 = sols.c1_6;
c2_6 = sols.c2_6;
c3_6 = sols.c3_6;
c4_6 = sols.c4_6;
c5_6 = sols.c5_6;
c6_6 = sols.c6_6;
c7_6 = sols.c7_6;

c1_7 = sols.c1_7;
c2_7 = sols.c2_7;
c3_7 = sols.c3_7;
c4_7 = sols.c4_7;
c5_7 = sols.c5_7;
c6_7 = sols.c6_7;
c7_7 = sols.c7_7;

save('wCoef_7state_condensed.mat','sols',...
    'c1_1', 'c1_2', 'c1_3', 'c1_4', 'c1_5', 'c1_6', 'c1_7', ...
    'c2_1', 'c2_2', 'c2_3', 'c2_4', 'c2_5', 'c2_6', 'c2_7',...
    'c3_1', 'c3_2', 'c3_3', 'c3_4', 'c3_5', 'c3_6', 'c3_7',...
    'c4_1', 'c4_2', 'c4_3', 'c4_4', 'c4_5', 'c4_6', 'c4_7',...
    'c5_1', 'c5_2', 'c5_3', 'c5_4', 'c5_5', 'c5_6', 'c5_7',...
    'c6_1', 'c6_2', 'c6_3', 'c6_4', 'c6_5', 'c6_6', 'c6_7',...
    'c7_1', 'c7_2', 'c7_3', 'c7_4', 'c7_5', 'c7_6', 'c7_7',...
    'P1_eq', 'P2_eq','P3_eq','P4_eq','P5_eq','P6_eq','P7_eq');

toc
%%
vars = sym.empty();
for j = 1:length(c)
    for i = 1:length(c)
        temp(i) = c(i,j);
        vars_temp(i) = temp(i);
    end
    vars = [vars, vars_temp]; % Builds column by column from c.
end



% This will not work, need to substitute values in by hand.
% for k = 1:length(vars)
%     vars(k) = sols.vars(k)
% end

