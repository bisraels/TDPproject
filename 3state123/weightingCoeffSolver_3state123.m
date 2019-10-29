
syms a b c %Components of the first eigenvector
syms x y z %Components of the second eigenvector
syms m_1 m_2 m_3 %Components of the first coefficient
syms n_1 n_2 n_3 %Components of the second coefficient
syms p1_eq p2_eq p3_eq %Equilibrium Populations

syms p_1 p_2 p_3  % coefficients for population starting in P1
syms f g h        % compoents of zeroth eigenvector

% % P i --> 3                          % Pji
% eqn1 = 1 == p3_eq + m_3*c + n_3*z;   % P33  
% eqn2 = 0 == p3_eq + m_2*c + n_2*z;   % P32
% eqn3 = 0 == p3_eq + m_1*c + n_1*z;   % P31
% 
% % P i --> 2
% eqn4 = 0 == p2_eq + m_3*b + n_3*y;   % P23
% eqn5 = 1 == p2_eq + m_2*b + n_2*y;   % P22
% eqn6 = 0 == p2_eq + m_1*b + n_1*y;   % P21
% 
% % P i --> 1
% eqn7 = 0 == p1_eq + m_3*a + n_3*x;   % P13
% eqn8 = 0 == p1_eq + m_2*a + n_2*x;   % P12
% eqn9 = 1 == p1_eq + m_1*a + n_1*x;   % P11


% The computer sees this an an inhomogenous differential equation which is
% hard to solve. If we make it easier by making it compatible with matrix
% form (break up the p_eq terms) it can solve it.


% Define p1, p2, p3 (another weighting coefficient) and f g h (components
% of 0th eigenvector)

% P i --> 3                          % Pji
eqn1 = 1 == p_3*h + m_3*c + n_3*z;   % P33  
eqn2 = 0 == p_2*h + m_2*c + n_2*z;   % P32
eqn3 = 0 == p_1*h + m_1*c + n_1*z;   % P31

% P i --> 2
eqn4 = 0 == p_3*g + m_3*b + n_3*y;   % P23
eqn5 = 1 == p_2*g + m_2*b + n_2*y;   % P22
eqn6 = 0 == p_1*g + m_1*b + n_1*y;   % P21

% P i --> 1
eqn7 = 0 == p_3*f + m_3*a + n_3*x;   % P13
eqn8 = 0 == p_2*f + m_2*a + n_2*x;   % P12
eqn9 = 1 == p_1*f + m_1*a + n_1*x;   % P11

p1_eq = p_1*f;
p2_eq = p_2*g;
p3_eq = p_3*h;

% eqn10 = p_1*f + p_2*g + p_3*h == 1;
eqn10 = p1_eq + p2_eq + p3_eq == 1;

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
vars = [m_1, m_2, m_3, n_1, n_2, n_3, p_1, p_2, p_3];
sol = solve(eqns,vars)

p1_eq = (sol.p_1)*f;
p2_eq = (sol.p_2)*g;
p3_eq = (sol.p_3)*h;

