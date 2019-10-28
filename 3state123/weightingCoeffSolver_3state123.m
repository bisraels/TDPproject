
syms a b c %Components of the first eigenvector
syms x y z %Components of the second eigenvector
syms m_1 m_2 m_3 %Components of the first coefficient
syms n_1 n_2 n_3 %Components of the second coefficient
syms p1_eq p2_eq p3_eq %Equilibrium Populations

eqn1 = 1 == p3_eq + m_3*c + n_3*z;
eqn2 = 0 == p3_eq + m_2*c + n_2*z;
eqn3 = 0 == p3_eq + m_1*c + n_1*z;

eqn4 = 0 == p2_eq + m_3*b + n_3*y;
eqn5 = 1 == p2_eq + m_2*b + n_2*y;
eqn6 = 0 == p2_eq + m_1*b + n_1*y;

eqn7 = 0 == p1_eq + m_3*a + n_3*x;
eqn8 = 0 == p1_eq + m_2*a + n_2*x;
eqn9 = 1 == p1_eq + m_1*a + n_1*x;

eqn10 = p1_eq + p2_eq + p3_eq == 1;

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
vars = [m_1, m_2, m_3, n_1, n_2, n_3];
sol = solve(eqns,vars)
