function [C4,C2_product,C4_diff] = FourPtTCF_3state_V2(tau1range,tau2,tau3range,A1,A2,A3,k12,k21,k23,k32)


% Initialize 4-pt. TCF matrices
C4_diff = zeros(length(tau1range),length(tau3range));
C4 = zeros(length(tau1range),length(tau3range));
C2_product = zeros(length(tau1range),length(tau3range));

% Define eigenvalues lam1 and lam2
c1 = 0.5*(k12+k21+k23+k32);
c2 = 0.5*sqrt(k12^2 + 2*k12*k21 - 2*k12*k23 - 2*k12*k32 + k21^2 + 2*k21*k23 - 2*k21*k32 + k23^2 + 2*k23*k32 + k32^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (c1 + c2)/k21 - (k12 + k21)/k21;
b = k12/k21 - (c1 + c2)/k21;
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (c1 - c2)/k21 - (k12 + k21)/k21;
y = k12/k21 - (c1 - c2)/k21;
z = 1;

% Define equilibrium values of populations
p1_eq = 1/(k21/k12 + 1 + k23/k32);
p0_eq = k21/k12*p1_eq;
p2_eq = k23/k32*p1_eq;


% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_2 = (-p1_eq - b/a*(1 - p2_eq))/(y - b/a*x);
n_1 = (1 - p1_eq + b/a*p2_eq)/(y - b/a*x);
n_0 = (-p1_eq + b/a*p2_eq)/(y - b/a*x);
m_2 = (1 - p2_eq - n_2*x)/a;
m_1 = (-p2_eq - n_1*x)/a;
m_0 = (-p2_eq - n_0*x)/a;

    t1 = tau1range;
    t2 = tau2;
    t3 = tau3range';

% Calculate TCF
path1 = p2_eq*A3*p2(2,t1)*A3*p2(2,t2)*A3.*p2(2,t3)*A3;
path2 = p2_eq*A3*p2(2,t1)*A3*p2(2,t2)*A3.*p1(2,t3)*A2;
path3 = p2_eq*A3*p2(2,t1)*A3*p2(2,t2)*A3.*p0(2,t3)*A1;
path4 = p2_eq*A3*p2(2,t1)*A3*p1(2,t2)*A2.*p2(1,t3)*A3;
path5 = p2_eq*A3*p2(2,t1)*A3*p1(2,t2)*A2.*p1(1,t3)*A2;
path6 = p2_eq*A3*p2(2,t1)*A3*p1(2,t2)*A2.*p0(1,t3)*A1;
path7 = p2_eq*A3*p2(2,t1)*A3*p0(2,t2)*A1.*p2(0,t3)*A3;
path8 = p2_eq*A3*p2(2,t1)*A3*p0(2,t2)*A1.*p1(0,t3)*A2;
path9 = p2_eq*A3*p2(2,t1)*A3*p0(2,t2)*A1.*p0(0,t3)*A1;

path10 = p2_eq*A3*p1(2,t1)*A2*p2(1,t2)*A3.*p2(2,t3)*A3;
path11 = p2_eq*A3*p1(2,t1)*A2*p2(1,t2)*A3.*p1(2,t3)*A2;
path12 = p2_eq*A3*p1(2,t1)*A2*p2(1,t2)*A3.*p0(2,t3)*A1;
path13 = p2_eq*A3*p1(2,t1)*A2*p1(1,t2)*A2.*p2(1,t3)*A3;
path14 = p2_eq*A3*p1(2,t1)*A2*p1(1,t2)*A2.*p1(1,t3)*A2;
path15 = p2_eq*A3*p1(2,t1)*A2*p1(1,t2)*A2.*p0(1,t3)*A1;
path16 = p2_eq*A3*p1(2,t1)*A2*p0(1,t2)*A1.*p2(0,t3)*A3;
path17 = p2_eq*A3*p1(2,t1)*A2*p0(1,t2)*A1.*p1(0,t3)*A2;
path18 = p2_eq*A3*p1(2,t1)*A2*p0(1,t2)*A1.*p0(0,t3)*A1;

path19 = p2_eq*A3*p0(2,t1)*A1*p2(0,t2)*A3.*p2(2,t3)*A3;
path20 = p2_eq*A3*p0(2,t1)*A1*p2(0,t2)*A3.*p1(2,t3)*A2;
path21 = p2_eq*A3*p0(2,t1)*A1*p2(0,t2)*A3.*p0(2,t3)*A1;
path22 = p2_eq*A3*p0(2,t1)*A1*p1(0,t2)*A2.*p2(1,t3)*A3;
path23 = p2_eq*A3*p0(2,t1)*A1*p1(0,t2)*A2.*p1(1,t3)*A2;
path24 = p2_eq*A3*p0(2,t1)*A1*p1(0,t2)*A2.*p0(1,t3)*A1;
path25 = p2_eq*A3*p0(2,t1)*A1*p0(0,t2)*A1.*p2(0,t3)*A3;
path26 = p2_eq*A3*p0(2,t1)*A1*p0(0,t2)*A1.*p1(0,t3)*A2;
path27 = p2_eq*A3*p0(2,t1)*A1*p0(0,t2)*A1.*p0(0,t3)*A1;

path28 = p1_eq*A2*p2(1,t1)*A3*p2(2,t2)*A3.*p2(2,t3)*A3;
path29 = p1_eq*A2*p2(1,t1)*A3*p2(2,t2)*A3.*p1(2,t3)*A2;
path30 = p1_eq*A2*p2(1,t1)*A3*p2(2,t2)*A3.*p0(2,t3)*A1;
path31 = p1_eq*A2*p2(1,t1)*A3*p1(2,t2)*A2.*p2(1,t3)*A3;
path32 = p1_eq*A2*p2(1,t1)*A3*p1(2,t2)*A2.*p1(1,t3)*A2;
path33 = p1_eq*A2*p2(1,t1)*A3*p1(2,t2)*A2.*p0(1,t3)*A1;
path34 = p1_eq*A2*p2(1,t1)*A3*p0(2,t2)*A1.*p2(0,t3)*A3;
path35 = p1_eq*A2*p2(1,t1)*A3*p0(2,t2)*A1.*p1(0,t3)*A2;
path36 = p1_eq*A2*p2(1,t1)*A3*p0(2,t2)*A1.*p0(0,t3)*A1;

path37 = p1_eq*A2*p1(1,t1)*A2*p2(1,t2)*A3.*p2(2,t3)*A3;
path38 = p1_eq*A2*p1(1,t1)*A2*p2(1,t2)*A3.*p1(2,t3)*A2;
path39 = p1_eq*A2*p1(1,t1)*A2*p2(1,t2)*A3.*p0(2,t3)*A1;
path40 = p1_eq*A2*p1(1,t1)*A2*p1(1,t2)*A2.*p2(1,t3)*A3;
path41 = p1_eq*A2*p1(1,t1)*A2*p1(1,t2)*A2.*p1(1,t3)*A2;
path42 = p1_eq*A2*p1(1,t1)*A2*p1(1,t2)*A2.*p0(1,t3)*A1;
path43 = p1_eq*A2*p1(1,t1)*A2*p0(1,t2)*A1.*p2(0,t3)*A3;
path44 = p1_eq*A2*p1(1,t1)*A2*p0(1,t2)*A1.*p1(0,t3)*A2;
path45 = p1_eq*A2*p1(1,t1)*A2*p0(1,t2)*A1.*p0(0,t3)*A1;

path46 = p1_eq*A2*p0(1,t1)*A1*p2(0,t2)*A3.*p2(2,t3)*A3;
path47 = p1_eq*A2*p0(1,t1)*A1*p2(0,t2)*A3.*p1(2,t3)*A2;
path48 = p1_eq*A2*p0(1,t1)*A1*p2(0,t2)*A3.*p0(2,t3)*A1;
path49 = p1_eq*A2*p0(1,t1)*A1*p1(0,t2)*A2.*p2(1,t3)*A3;
path50 = p1_eq*A2*p0(1,t1)*A1*p1(0,t2)*A2.*p1(1,t3)*A2;
path51 = p1_eq*A2*p0(1,t1)*A1*p1(0,t2)*A2.*p0(1,t3)*A1;
path52 = p1_eq*A2*p0(1,t1)*A1*p0(0,t2)*A1.*p2(0,t3)*A3;
path53 = p1_eq*A2*p0(1,t1)*A1*p0(0,t2)*A1.*p1(0,t3)*A2;
path54 = p1_eq*A2*p0(1,t1)*A1*p0(0,t2)*A1.*p0(0,t3)*A1;

path55 = p0_eq*A1*p2(0,t1)*A3*p2(2,t2)*A3.*p2(2,t3)*A3;
path56 = p0_eq*A1*p2(0,t1)*A3*p2(2,t2)*A3.*p1(2,t3)*A2;
path57 = p0_eq*A1*p2(0,t1)*A3*p2(2,t2)*A3.*p0(2,t3)*A1;
path58 = p0_eq*A1*p2(0,t1)*A3*p1(2,t2)*A2.*p2(1,t3)*A3;
path59 = p0_eq*A1*p2(0,t1)*A3*p1(2,t2)*A2.*p1(1,t3)*A2;
path60 = p0_eq*A1*p2(0,t1)*A3*p1(2,t2)*A2.*p0(1,t3)*A1;
path61 = p0_eq*A1*p2(0,t1)*A3*p0(2,t2)*A1.*p2(0,t3)*A3;
path62 = p0_eq*A1*p2(0,t1)*A3*p0(2,t2)*A1.*p1(0,t3)*A2;
path63 = p0_eq*A1*p2(0,t1)*A3*p0(2,t2)*A1.*p0(0,t3)*A1;

path64 = p0_eq*A1*p1(0,t1)*A2*p2(1,t2)*A3.*p2(2,t3)*A3;
path65 = p0_eq*A1*p1(0,t1)*A2*p2(1,t2)*A3.*p1(2,t3)*A2;
path66 = p0_eq*A1*p1(0,t1)*A2*p2(1,t2)*A3.*p0(2,t3)*A1;
path67 = p0_eq*A1*p1(0,t1)*A2*p1(1,t2)*A2.*p2(1,t3)*A3;
path68 = p0_eq*A1*p1(0,t1)*A2*p1(1,t2)*A2.*p1(1,t3)*A2;
path69 = p0_eq*A1*p1(0,t1)*A2*p1(1,t2)*A2.*p0(1,t3)*A1;
path70 = p0_eq*A1*p1(0,t1)*A2*p0(1,t2)*A1.*p2(0,t3)*A3;
path71 = p0_eq*A1*p1(0,t1)*A2*p0(1,t2)*A1.*p1(0,t3)*A2;
path72 = p0_eq*A1*p1(0,t1)*A2*p0(1,t2)*A1.*p0(0,t3)*A1;

path73 = p0_eq*A1*p0(0,t1)*A1*p2(0,t2)*A3.*p2(2,t3)*A3;
path74 = p0_eq*A1*p0(0,t1)*A1*p2(0,t2)*A3.*p1(2,t3)*A2;
path75 = p0_eq*A1*p0(0,t1)*A1*p2(0,t2)*A3.*p0(2,t3)*A1;
path76 = p0_eq*A1*p0(0,t1)*A1*p1(0,t2)*A2.*p2(1,t3)*A3;
path77 = p0_eq*A1*p0(0,t1)*A1*p1(0,t2)*A2.*p1(1,t3)*A2;
path78 = p0_eq*A1*p0(0,t1)*A1*p1(0,t2)*A2.*p0(1,t3)*A1;
path79 = p0_eq*A1*p0(0,t1)*A1*p0(0,t2)*A1.*p2(0,t3)*A3;
path80 = p0_eq*A1*p0(0,t1)*A1*p0(0,t2)*A1.*p1(0,t3)*A2;
path81 = p0_eq*A1*p0(0,t1)*A1*p0(0,t2)*A1.*p0(0,t3)*A1;

kappa2 = path1 + path2 + path3 + path4 + path5 + path6 + path7 + path8 + path9 + path10 + path11 + path12 + path13 + path14 + path15 + path16 + path17 + path18 + path19 + path20;
kappa2 = kappa2 + path21 + path22 + path23 + path24 + path25 + path26 + path27 + path28 + path29 + path30 + path31 + path32 + path33 + path34 + path35 + path36 + path37 + path38 + path39 + path40;
kappa2 = kappa2 + path41 + path42 + path43 + path44 + path45 + path46 + path47 + path48 + path49 + path50 + path51 + path52 + path53 + path54 + path55 + path56 + path57 + path58 + path59 + path60;
C4 = kappa2 + path61 + path62 + path63 + path64 + path65 + path66 + path67 + path68 + path69 + path70 + path71 + path72 + path73 + path74 + path75 + path76 + path77 + path78 + path79 + path80 + path81;
theta_1 = (p2_eq*A3*(A3*p2(2,t1) + A2*p1(2,t1) + A1*p0(2,t1)) + p1_eq*A2*(A3*p2(1,t1) + A2*p1(1,t1) + A1*p0(1,t1)) + p0_eq*A1*(A3*p2(0,t1) + A2*p1(0,t1) + A1*p0(0,t1)));
theta_2 = (p2_eq*A3*(A3*p2(2,t3) + A2*p1(2,t3) + A1*p0(2,t3)) + p1_eq*A2*(A3*p2(1,t3) + A2*p1(1,t3) + A1*p0(1,t3)) + p0_eq*A1*(A3*p2(0,t3) + A2*p1(0,t3) + A1*p0(0,t3)));
C2_product = theta_1.*theta_2';

C4_diff = C4 - C2_product;

function p0 = p0(init,time)
    if init == 2
        m = m_2;
        n = n_2;
    elseif init == 1
        m = m_1;
        n = n_1;
    else
        m = m_0;
        n = n_0;
    end

    p0 = p0_eq + m*c*exp(-lam1*time) + n*z*exp(-lam2*time);
end

function p1 = p1(init,time)
    if init == 2
        m = m_2;
        n = n_2;
    elseif init == 1
        m = m_1;
        n = n_1;
    else
        m = m_0;
        n = n_0;
    end

    p1 = p1_eq + m*b*exp(-lam1*time) + n*y*exp(-lam2*time);
end

function p2 = p2(init,time)
    if init == 2
        m = m_2;
        n = n_2;
    elseif init == 1
        m = m_1;
        n = n_1;
    else
        m = m_0;
        n = n_0;
    end

    p2 = p2_eq + m*a*exp(-lam1*time) + n*x*exp(-lam2*time);
end

end