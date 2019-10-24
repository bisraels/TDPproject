syms k12 k21 k23 k32 k13 k31
syms P1eq P2eq P3eq

% Define the Rate Matrix K
K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];

P = [P1eq; P2eq; P3eq];

dP_t = [0; 0 ; 0];
%@Equilibrium: dP(t)/dt = K*P(t)
eqn = dP_t == K*P;

sol = solve([eqn(1),eqn(2),eqn(3)],[P1eq,P2eq,P3eq]);

P1eqSol = sol.P1eq;
P2eqSol = sol.P2eq;
P3eqSol = sol.P3eq;