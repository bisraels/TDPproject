syms k12 k21 k23 k32 k13 k31
syms P1eq P2eq P3eq

% Define the Rate Matrix K
K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];

P = [P1eq; P2eq; P3eq];

dP_t = [0; 0 ; 0];
%@Equilibrium: dP(t)/dt = K*P(t)
eqn = K*P == dP_t

eqn4 = P1eq + P2eq + P3eq == 1;
eqns = [eqn(1),eqn(2),eqn(3),eqn4];
% % % OUTPUT::::
% % % eqn =
% % %  
% % %  P2eq*k21 + P3eq*k31 - P1eq*(k12 + k13) == 0
% % %  P1eq*k12 + P3eq*k32 - P2eq*(k21 + k23) == 0
% % %  P1eq*k13 + P2eq*k23 - P3eq*(k31 + k32) == 0

% [A,B] = equationsToMatrix([eqn(1),eqn(2),eqn(3)], [P1eq,P2eq,P3eq]);
% X = linsolve(A,B)


sol = solve(eqns,[P1eq,P2eq,P3eq]);
P1eqSol = sol.P1eq
P2eqSol = sol.P2eq
P3eqSol = sol.P3eq

% %Solve the sysyem of equations for one of three variables
% P2eq_sol = solve(eqn(1),P2eq)
% %Plug the solution of the equation back into the overall equation
% eqn2 = subs(eqn,P2eq,P2eq_sol)
% %solve the system of equations for another variable
% P3eq_sol = solve(eqn2(2),P3eq)
% eqn3 = subs(eqn2,P3eq,P3eq_sol)
% P1eq_sol = solve(eqn3(1),P1eq)
%  