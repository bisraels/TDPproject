
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
eq10 = p1_eq + p2_eq + p3_eq == 1;

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
vars = [m_1, m_2, m_3, n_1, n_2, n_3, p_1, p_2, p_3];
sol = solve(eqns,vars)

p1_eq = (sol.p_1)*f;
p2_eq = (sol.p_2)*g;
p3_eq = (sol.p_3)*h;


%%


% eqn9 = 1 == h0 * d + m_1 * a + n_1 * x;   % P11
% eqn6 = 0 == h0 * e + m_1 * b + n_1 * y;   % P21
% eqn3 = 0 == h0 * f + m_1 * c + n_1 * z;   % P31
% 
% sols = solve([eqn9, eqn6, eqn3],[m_1, n_1, h0])

%% Try for 7 state

syms c1_1 c1_2 c1_3 c1_4 c1_5 c1_6 c1_7
syms c2_1 c2_2 c2_3 c2_4 c2_5 c2_6 c2_7
syms c3_1 c3_2 c3_3 c3_4 c3_5 c3_6 c3_7
syms c4_1 c4_2 c4_3 c4_4 c4_5 c4_6 c4_7
syms c5_1 c5_2 c5_3 c5_4 c5_5 c5_6 c5_7
syms c6_1 c6_2 c6_3 c6_4 c6_5 c6_6 c6_7
syms c7_1 c7_2 c7_3 c7_4 c7_5 c7_6 c7_7

syms v1_1 v1_2 v1_3 v1_4 v1_5 v1_6 v1_7
syms v2_1 v2_2 v2_3 v2_4 v2_5 v2_6 v2_7
syms v3_1 v3_2 v3_3 v3_4 v3_5 v3_6 v3_7
syms v4_1 v4_2 v4_3 v4_4 v4_5 v4_6 v4_7
syms v5_1 v5_2 v5_3 v5_4 v5_5 v5_6 v5_7
syms v6_1 v6_2 v6_3 v6_4 v6_5 v6_6 v6_7
syms v7_1 v7_2 v7_3 v7_4 v7_5 v7_6 v7_7


% Coefficients: cm_n with  m = vector, n = condition
% Vectors: vm_n with m = vector, n = component

% Note: c1_n * v1_n = pn_eq

  
% Pji :  P i --> j

% Pj1: condition = 1
% P11
eqn1_1 = 1 == (c1_1 * v1_1) + (c2_1 * v2_1) + (c3_1 * v3_1) + (c4_1 * v4_1) + (c5_1 * v5_1) + (c6_1 * v6_1) + (c7_1 * v7_1);  
% P21
eqn1_2 = 0 == (c1_1 * v1_2) + (c2_1 * v2_2) + (c3_1 * v3_2) + (c4_1 * v4_2) + (c5_1 * v5_2) + (c6_1 * v6_2) + (c7_1 * v7_2);
% P31
eqn1_3 = 0 == (c1_1 * v1_3) + (c2_1 * v2_3) + (c3_1 * v3_3) + (c4_1 * v4_3) + (c5_1 * v5_3) + (c6_1 * v6_3) + (c7_1 * v7_3);
% P41
eqn1_4 = 0 == (c1_1 * v1_4) + (c2_1 * v2_4) + (c3_1 * v3_4) + (c4_1 * v4_4) + (c5_1 * v5_4) + (c6_1 * v6_4) + (c7_1 * v7_4);
% P51
eqn1_5 = 0 == (c1_1 * v1_5) + (c2_1 * v2_5) + (c3_1 * v3_5) + (c4_1 * v4_5) + (c5_1 * v5_5) + (c6_1 * v6_5) + (c7_1 * v7_5);
% P61
eqn1_6 = 0 == (c1_1 * v1_6) + (c2_1 * v2_6) + (c3_1 * v3_6) + (c4_1 * v4_6) + (c5_1 * v5_6) + (c6_1 * v6_6) + (c7_1 * v7_6);
% P71
eqn1_7 = 0 == (c1_1 * v1_7) + (c2_1 * v2_7) + (c3_1 * v3_7) + (c4_1 * v4_7) + (c5_1 * v5_7) + (c6_1 * v6_7) + (c7_1 * v7_7);


% Pj2: condtion = 2
% P12
eqn2_1 = 0 == (c1_2 * v1_1) + (c2_2 * v2_1) + (c3_2 * v3_1) + (c4_2 * v4_1) + (c5_2 * v5_1) + (c6_2 * v6_1) + (c7_2 * v7_1);  
% P22
eqn2_2 = 1 == (c1_2 * v1_2) + (c2_2 * v2_2) + (c3_2 * v3_2) + (c4_2 * v4_2) + (c5_2 * v5_2) + (c6_2 * v6_2) + (c7_2 * v7_2);
% P32
eqn2_3 = 0 == (c1_2 * v1_3) + (c2_2 * v2_3) + (c3_2 * v3_3) + (c4_2 * v4_3) + (c5_2 * v5_3) + (c6_2 * v6_3) + (c7_2 * v7_3);
% P42
eqn2_4 = 0 == (c1_2 * v1_4) + (c2_2 * v2_4) + (c3_2 * v3_4) + (c4_2 * v4_4) + (c5_2 * v5_4) + (c6_2 * v6_4) + (c7_2 * v7_4);
% P52
eqn2_5 = 0 == (c1_2 * v1_5) + (c2_2 * v2_5) + (c3_2 * v3_5) + (c4_2 * v4_5) + (c5_2 * v5_5) + (c6_2 * v6_5) + (c7_2 * v7_5);
% P62
eqn2_6 = 0 == (c1_2 * v1_6) + (c2_2 * v2_6) + (c3_2 * v3_6) + (c4_2 * v4_6) + (c5_2 * v5_6) + (c6_2 * v6_6) + (c7_2 * v7_6);
% P72
eqn2_7 = 0 == (c1_2 * v1_7) + (c2_2 * v2_7) + (c3_2 * v3_7) + (c4_2 * v4_7) + (c5_2 * v5_7) + (c6_2 * v6_7) + (c7_2 * v7_7);


% Pj3: condtion = 3
% P13
eqn3_1 = 0 == (c1_3 * v1_1) + (c2_3 * v2_1) + (c3_3 * v3_1) + (c4_3 * v4_1) + (c5_3 * v5_1) + (c6_3 * v6_1) + (c7_3 * v7_1);  
% P23
eqn3_2 = 0 == (c1_3 * v1_2) + (c2_3 * v2_2) + (c3_3 * v3_2) + (c4_3 * v4_2) + (c5_3 * v5_2) + (c6_3 * v6_2) + (c7_3 * v7_2);
% P33
eqn3_3 = 1 == (c1_3 * v1_3) + (c2_3 * v2_3) + (c3_3 * v3_3) + (c4_3 * v4_3) + (c5_3 * v5_3) + (c6_3 * v6_3) + (c7_3 * v7_3);
% P43
eqn3_4 = 0 == (c1_3 * v1_4) + (c2_3 * v2_4) + (c3_3 * v3_4) + (c4_3 * v4_4) + (c5_3 * v5_4) + (c6_3 * v6_4) + (c7_3 * v7_4);
% P53
eqn3_5 = 0 == (c1_3 * v1_5) + (c2_3 * v2_5) + (c3_3 * v3_5) + (c4_3 * v4_5) + (c5_3 * v5_5) + (c6_3 * v6_5) + (c7_3 * v7_5);
% P63
eqn3_6 = 0 == (c1_3 * v1_6) + (c2_3 * v2_6) + (c3_3 * v3_6) + (c4_3 * v4_6) + (c5_3 * v5_6) + (c6_3 * v6_6) + (c7_3 * v7_6);
% P73
eqn3_7 = 0 == (c1_3 * v1_7) + (c2_3 * v2_7) + (c3_3 * v3_7) + (c4_3 * v4_7) + (c5_3 * v5_7) + (c6_3 * v6_7) + (c7_3 * v7_7);


% Pj4: condtion = 4
% P14
eqn4_1 = 0 == (c1_4 * v1_1) + (c2_4 * v2_1) + (c3_4 * v3_1) + (c4_4 * v4_1) + (c5_4 * v5_1) + (c6_4 * v6_1) + (c7_4 * v7_1);  
% P24
eqn4_2 = 0 == (c1_4 * v1_2) + (c2_4 * v2_2) + (c3_4 * v3_2) + (c4_4 * v4_2) + (c5_4 * v5_2) + (c6_4 * v6_2) + (c7_4 * v7_2);
% P34
eqn4_3 = 0 == (c1_4 * v1_3) + (c2_4 * v2_3) + (c3_4 * v3_3) + (c4_4 * v4_3) + (c5_4 * v5_3) + (c6_4 * v6_3) + (c7_4 * v7_3);
% P44
eqn4_4 = 1 == (c1_4 * v1_4) + (c2_4 * v2_4) + (c3_4 * v3_4) + (c4_4 * v4_4) + (c5_4 * v5_4) + (c6_4 * v6_4) + (c7_4 * v7_4);
% P54
eqn4_5 = 0 == (c1_4 * v1_5) + (c2_4 * v2_5) + (c3_4 * v3_5) + (c4_4 * v4_5) + (c5_4 * v5_5) + (c6_4 * v6_5) + (c7_4 * v7_5);
% P64
eqn4_6 = 0 == (c1_4 * v1_6) + (c2_4 * v2_6) + (c3_4 * v3_6) + (c4_4 * v4_6) + (c5_4 * v5_6) + (c6_4 * v6_6) + (c7_4 * v7_6);
% P74
eqn4_7 = 0 == (c1_4 * v1_7) + (c2_4 * v2_7) + (c3_4 * v3_7) + (c4_4 * v4_7) + (c5_4 * v5_7) + (c6_4 * v6_7) + (c7_4 * v7_7);

% Pj5: condtion = 5
% P15
eqn5_1 = 0 == (c1_5 * v1_1) + (c2_5 * v2_1) + (c3_5 * v3_1) + (c4_5 * v4_1) + (c5_5 * v5_1) + (c6_5 * v6_1) + (c7_5 * v7_1);  
% P25
eqn5_2 = 0 == (c1_5 * v1_2) + (c2_5 * v2_2) + (c3_5 * v3_2) + (c4_5 * v4_2) + (c5_5 * v5_2) + (c6_5 * v6_2) + (c7_5 * v7_2);
% P35
eqn5_3 = 0 == (c1_5 * v1_3) + (c2_5 * v2_3) + (c3_5 * v3_3) + (c4_5 * v4_3) + (c5_5 * v5_3) + (c6_5 * v6_3) + (c7_5 * v7_3);
% P45
eqn5_4 = 0 == (c1_5 * v1_4) + (c2_5 * v2_4) + (c3_5 * v3_4) + (c4_5 * v4_4) + (c5_5 * v5_4) + (c6_5 * v6_4) + (c7_5 * v7_4);
% P55
eqn5_5 = 1 == (c1_5 * v1_5) + (c2_5 * v2_5) + (c3_5 * v3_5) + (c4_5 * v4_5) + (c5_5 * v5_5) + (c6_5 * v6_5) + (c7_5 * v7_5);
% P65
eqn5_6 = 0 == (c1_5 * v1_6) + (c2_5 * v2_6) + (c3_5 * v3_6) + (c4_5 * v4_6) + (c5_5 * v5_6) + (c6_5 * v6_6) + (c7_5 * v7_6);
% P75
eqn5_7 = 0 == (c1_5 * v1_7) + (c2_5 * v2_7) + (c3_5 * v3_7) + (c4_5 * v4_7) + (c5_5 * v5_7) + (c6_5 * v6_7) + (c7_5 * v7_7);

% Pj6: condtion = 6
% P16
eqn6_1 = 0 == (c1_6 * v1_1) + (c2_6 * v2_1) + (c3_6 * v3_1) + (c4_6 * v4_1) + (c5_6 * v5_1) + (c6_6 * v6_1) + (c7_6 * v7_1);  
% P26
eqn6_2 = 0 == (c1_6 * v1_2) + (c2_6 * v2_2) + (c3_6 * v3_2) + (c4_6 * v4_2) + (c5_6 * v5_2) + (c6_6 * v6_2) + (c7_6 * v7_2);
% P36
eqn6_3 = 0 == (c1_6 * v1_3) + (c2_6 * v2_3) + (c3_6 * v3_3) + (c4_6 * v4_3) + (c5_6 * v5_3) + (c6_6 * v6_3) + (c7_6 * v7_3);
% P46
eqn6_4 = 0 == (c1_6 * v1_4) + (c2_6 * v2_4) + (c3_6 * v3_4) + (c4_6 * v4_4) + (c5_6 * v5_4) + (c6_6 * v6_4) + (c7_6 * v7_4);
% P56
eqn6_5 = 0 == (c1_6 * v1_5) + (c2_6 * v2_5) + (c3_6 * v3_5) + (c4_6 * v4_5) + (c5_6 * v5_5) + (c6_6 * v6_5) + (c7_6 * v7_5);
% P66
eqn6_6 = 1 == (c1_6 * v1_6) + (c2_6 * v2_6) + (c3_6 * v3_6) + (c4_6 * v4_6) + (c5_6 * v5_6) + (c6_6 * v6_6) + (c7_6 * v7_6);
% P76
eqn6_7 = 0 == (c1_6 * v1_7) + (c2_6 * v2_7) + (c3_6 * v3_7) + (c4_6 * v4_7) + (c5_6 * v5_7) + (c6_6 * v6_7) + (c7_6 * v7_7);

% Pj7: condtion = 7
% P17
eqn7_1 = 0 == (c1_7 * v1_1) + (c2_7 * v2_1) + (c3_7 * v3_1) + (c4_7 * v4_1) + (c5_7 * v5_1) + (c6_7 * v6_1) + (c7_6 * v7_1);  
% P27
eqn7_2 = 0 == (c1_7 * v1_2) + (c2_7 * v2_2) + (c3_7 * v3_2) + (c4_7 * v4_2) + (c5_7 * v5_2) + (c6_7 * v6_2) + (c7_6 * v7_2);
% P37
eqn7_3 = 0 == (c1_7 * v1_3) + (c2_7 * v2_3) + (c3_7 * v3_3) + (c4_7 * v4_3) + (c5_7 * v5_3) + (c6_7 * v6_3) + (c7_6 * v7_3);
% P47
eqn7_4 = 0 == (c1_7 * v1_4) + (c2_7 * v2_4) + (c3_7 * v3_4) + (c4_7 * v4_4) + (c5_7 * v5_4) + (c6_7 * v6_4) + (c7_6 * v7_4);
% P57
eqn7_5 = 0 == (c1_7 * v1_5) + (c2_7 * v2_5) + (c3_7 * v3_5) + (c4_7 * v4_5) + (c5_7 * v5_5) + (c6_7 * v6_5) + (c7_6 * v7_5);
% P67
eqn7_6 = 0 == (c1_7 * v1_6) + (c2_7 * v2_6) + (c3_7 * v3_6) + (c4_7 * v4_6) + (c5_7 * v5_6) + (c6_7 * v6_6) + (c7_6 * v7_6);
% P77
eqn7_7 = 1 == (c1_7 * v1_7) + (c2_7 * v2_7) + (c3_7 * v3_7) + (c4_7 * v4_7) + (c5_7 * v5_7) + (c6_7 * v6_7) + (c7_6 * v7_7);


P1_eq = (c1_1 * v1_1); % from P11
P2_eq = (c1_2 * v1_2); % from P22
P3_eq = (c1_3 * v1_3); % from P33
P4_eq = (c1_4 * v1_4); % from P44
P5_eq = (c1_5 * v1_5); % from P55
P6_eq = (c1_6 * v1_6); % from P66
P7_eq = (c1_7 * v1_7); % from P77

eqn_sum = P1_eq + P2_eq + P3_eq + P4_eq + P5_eq + P6_eq + P7_eq == 1;


% Compact form to make sure they are all  there
% eqns = [eqn1_1, eqn1_2, eqn1_3, eqn1_4, eqn1_5, eqn1_6, eqn1_7,
%         eqn2_1, eqn2_2, eqn2_3, eqn2_4, eqn2_5, eqn2_6, eqn2_7,
%         eqn3_1, eqn3_2, eqn3_3, eqn3_4, eqn3_5, eqn3_6, eqn3_7,
%         eqn4_1, eqn4_2, eqn4_3, eqn4_4, eqn4_5, eqn4_6, eqn4_7,
%         eqn5_1, eqn5_2, eqn5_3, eqn5_4, eqn5_5, eqn5_6, eqn5_7,
%         eqn6_1, eqn6_2, eqn6_3, eqn6_4, eqn6_5, eqn6_6, eqn6_7,
%         eqn7_1, eqn7_2, eqn7_3, eqn7_4, eqn7_5, eqn7_6, eqn7_7,
%         eqn_sum];

% Expanded form  to make  computer happy.
eqns = [eqn1_1, eqn1_2, eqn1_3, eqn1_4, eqn1_5, eqn1_6, eqn1_7, eqn2_1, eqn2_2, eqn2_3, eqn2_4, eqn2_5, eqn2_6, eqn2_7, eqn3_1, eqn3_2, eqn3_3, eqn3_4, eqn3_5, eqn3_6, eqn3_7, eqn4_1, eqn4_2, eqn4_3, eqn4_4, eqn4_5, eqn4_6, eqn4_7, eqn5_1, eqn5_2, eqn5_3, eqn5_4, eqn5_5, eqn5_6, eqn5_7, eqn6_1, eqn6_2, eqn6_3, eqn6_4, eqn6_5, eqn6_6, eqn6_7, eqn7_1, eqn7_2, eqn7_3, eqn7_4, eqn7_5, eqn7_6, eqn7_7, eqn_sum];
        
% Compact form to make sure they are all  there
% vars = [c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7,
%         c2_1, c2_2, c2_3, c2_4, c2_5, c2_6, c2_7,
%         c3_1, c3_2, c3_3, c3_4, c3_5, c3_6, c3_7,
%         c4_1, c4_2, c4_3, c4_4, c4_5, c4_6, c4_7,
%         c5_1, c5_2, c5_3, c5_4, c5_5, c5_6, c5_7,
%         c6_1, c6_2, c6_3, c6_4, c6_5, c6_6, c6_7,
%         c7_1, c7_2, c7_3, c7_4, c7_5, c7_6, c7_7];

% Expanded form to make computer happy.
vars = [c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c2_1, c2_2, c2_3, c2_4, c2_5, c2_6, c2_7, c3_1, c3_2, c3_3, c3_4, c3_5, c3_6, c3_7, c4_1, c4_2, c4_3, c4_4, c4_5, c4_6, c4_7, c5_1, c5_2, c5_3, c5_4, c5_5, c5_6, c5_7, c6_1, c6_2, c6_3, c6_4, c6_5, c6_6, c6_7, c7_1, c7_2, c7_3, c7_4, c7_5, c7_6, c7_7];
      
sols = solve(eqns, vars)


% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * v1_1); 
P2_eq = (sols.c1_2 * v1_2); 
P3_eq = (sols.c1_3 * v1_3); 
P4_eq = (sols.c1_4 * v1_4); 
p5_eq = (sols.c1_5 * v1_5);
P6_eq = (sols.c1_6 * v1_6); 
P7_eq = (sols.c1_7 * v1_7); 