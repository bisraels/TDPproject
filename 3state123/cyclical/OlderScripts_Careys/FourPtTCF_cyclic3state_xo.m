function [p0_eq,p1_eq,p2_eq,lam1,lam2] = FourPtTCF_cyclic3state_xo(tau1range,tau2,A0,A1,A2,k01,k10,k12,k21,k20)

tau3range = tau1range;

k02 = k20*k12*k01/(k21*k10); % This relation forces detailed balance to be satisfied.

% Initialize 4-pt. TCF matrices
tcf = zeros(length(tau1range),length(tau3range));
kappa = zeros(length(tau1range),length(tau3range));
theta = zeros(length(tau1range),length(tau3range));

% Define eigenvalues lam1 and lam2
c1 = 0.5*(k01+k10+k12+k21+k02+k20);
c2 = 0.5*sqrt(k01^2 + 2*k01*k02 + 2*k01*k10 - 2*k01*k12 - 2*k01*k20 - 2*k01*k21 + k02^2 - 2*k02*k10 - 2*k02*k12 + 2*k02*k20 - k02*k21 + k10^2 + 2*k10*k12 - 2*k10*k20 - 2*k10*k21 + k12^2 - 2*k12*k20 + 2*k12*k21 + k20^2 + 2*k20*k21 + k21^2);
lam1 = c1 + c2;
lam2 = c1 - c2;

% Define eigenvector components corresponding to lam1: v1 = [a,b,c], in the
% basis [p2,p1,p0], where p# is the population in state #.
a = (k12 + k20 + k21 - c1 - c2)/(k02 - k12);
b = (c1 + c2 - k02 - k20 - k21)/(k02 - k12);
c = 1;

% Define eigenvector components corresponding to lam2: v2 = [x,y,z].
x = (k12 + k20 + k21 - c1 + c2)/(k02 - k12);
y = (c1 - c2 - k02 - k20 - k21)/(k02 - k12);
z = 1;

% Define equilibrium values of populations
p0_eq = 1/(1 + (k01 + k02)/k10 + (1 - k20/k10)*((k10 + k12)*(k01 + k02) - k01*k10)/((k10 + k12)*k20 + k21*k10));
p2_eq = ((k10 + k12)*(k01 + k02) - k01*k10)/((k10 + k12)*k20 + k21*k10)*p0_eq;
p1_eq = ((k01 + k02)*p0_eq - k20*p2_eq)/k10;

% Define constants needed for differential equation solutions, for various
% initial values of the population. m_2 is the constant "m" if the
% population of state 2 begins at 1 (others would then be zero).
n_2 = (-p1_eq + b/a*p0_eq)/(y - b/a*x);
n_1 = (1 - p1_eq + b/a*p0_eq)/(y - b/a*x);
n_0 = (-p1_eq - b/a*(1 - p0_eq))/(y - b/a*x);
m_2 = (-p0_eq - n_2*x)/a;
m_1 = (-p0_eq - n_1*x)/a;
m_0 = (1 - p0_eq - n_0*x)/a;

t1 = tau1range;
t2 = tau2;
t3 = tau3range';

% Subtract mean values
Amean = p0_eq*A0 + p1_eq*A1 + p2_eq*A2;
A0 = A0 - Amean;
A1 = A1 - Amean;
A2 = A2 - Amean;

% Calculate TCF
path1 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2*p2(2,t3)*A2;
path2 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2*p1(2,t3)*A1;
path3 = p2_eq*A2*p2(2,t1)*A2*p2(2,t2)*A2*p0(2,t3)*A0;
path4 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1*p2(1,t3)*A2;
path5 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1*p1(1,t3)*A1;
path6 = p2_eq*A2*p2(2,t1)*A2*p1(2,t2)*A1*p0(1,t3)*A0;
path7 = p2_eq*A2*p2(2,t1)*A2*p0(2,t2)*A0*p2(0,t3)*A2;
path8 = p2_eq*A2*p2(2,t1)*A2*p0(2,t2)*A0*p1(0,t3)*A1;
path9 = p2_eq*A2*p2(2,t1)*A2*p0(2,t2)*A0*p0(0,t3)*A0;

path10 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2*p2(2,t3)*A2;
path11 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2*p1(2,t3)*A1;
path12 = p2_eq*A2*p1(2,t1)*A1*p2(1,t2)*A2*p0(2,t3)*A0;
path13 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1*p2(1,t3)*A2;
path14 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1*p1(1,t3)*A1;
path15 = p2_eq*A2*p1(2,t1)*A1*p1(1,t2)*A1*p0(1,t3)*A0;
path16 = p2_eq*A2*p1(2,t1)*A1*p0(1,t2)*A0*p2(0,t3)*A2;
path17 = p2_eq*A2*p1(2,t1)*A1*p0(1,t2)*A0*p1(0,t3)*A1;
path18 = p2_eq*A2*p1(2,t1)*A1*p0(1,t2)*A0*p0(0,t3)*A0;

path19 = p2_eq*A2*p0(2,t1)*A0*p2(0,t2)*A2*p2(2,t3)*A2;
path20 = p2_eq*A2*p0(2,t1)*A0*p2(0,t2)*A2*p1(2,t3)*A1;
path21 = p2_eq*A2*p0(2,t1)*A0*p2(0,t2)*A2*p0(2,t3)*A0;
path22 = p2_eq*A2*p0(2,t1)*A0*p1(0,t2)*A1*p2(1,t3)*A2;
path23 = p2_eq*A2*p0(2,t1)*A0*p1(0,t2)*A1*p1(1,t3)*A1;
path24 = p2_eq*A2*p0(2,t1)*A0*p1(0,t2)*A1*p0(1,t3)*A0;
path25 = p2_eq*A2*p0(2,t1)*A0*p0(0,t2)*A0*p2(0,t3)*A2;
path26 = p2_eq*A2*p0(2,t1)*A0*p0(0,t2)*A0*p1(0,t3)*A1;
path27 = p2_eq*A2*p0(2,t1)*A0*p0(0,t2)*A0*p0(0,t3)*A0;

path28 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2*p2(2,t3)*A2;
path29 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2*p1(2,t3)*A1;
path30 = p1_eq*A1*p2(1,t1)*A2*p2(2,t2)*A2*p0(2,t3)*A0;
path31 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1*p2(1,t3)*A2;
path32 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1*p1(1,t3)*A1;
path33 = p1_eq*A1*p2(1,t1)*A2*p1(2,t2)*A1*p0(1,t3)*A0;
path34 = p1_eq*A1*p2(1,t1)*A2*p0(2,t2)*A0*p2(0,t3)*A2;
path35 = p1_eq*A1*p2(1,t1)*A2*p0(2,t2)*A0*p1(0,t3)*A1;
path36 = p1_eq*A1*p2(1,t1)*A2*p0(2,t2)*A0*p0(0,t3)*A0;

path37 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2*p2(2,t3)*A2;
path38 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2*p1(2,t3)*A1;
path39 = p1_eq*A1*p1(1,t1)*A1*p2(1,t2)*A2*p0(2,t3)*A0;
path40 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1*p2(1,t3)*A2;
path41 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1*p1(1,t3)*A1;
path42 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1*p0(1,t3)*A0;
path43 = p1_eq*A1*p1(1,t1)*A1*p0(1,t2)*A0*p2(0,t3)*A2;
path44 = p1_eq*A1*p1(1,t1)*A1*p0(1,t2)*A0*p1(0,t3)*A1;
path45 = p1_eq*A1*p1(1,t1)*A1*p0(1,t2)*A0*p0(0,t3)*A0;

path46 = p1_eq*A1*p0(1,t1)*A0*p2(0,t2)*A2*p2(2,t3)*A2;
path47 = p1_eq*A1*p0(1,t1)*A0*p2(0,t2)*A2*p1(2,t3)*A1;
path48 = p1_eq*A1*p0(1,t1)*A0*p2(0,t2)*A2*p0(2,t3)*A0;
path49 = p1_eq*A1*p0(1,t1)*A0*p1(0,t2)*A1*p2(1,t3)*A2;
path50 = p1_eq*A1*p0(1,t1)*A0*p1(0,t2)*A1*p1(1,t3)*A1;
path51 = p1_eq*A1*p0(1,t1)*A0*p1(0,t2)*A1*p0(1,t3)*A0;
path52 = p1_eq*A1*p0(1,t1)*A0*p0(0,t2)*A0*p2(0,t3)*A2;
path53 = p1_eq*A1*p0(1,t1)*A0*p0(0,t2)*A0*p1(0,t3)*A1;
path54 = p1_eq*A1*p0(1,t1)*A0*p0(0,t2)*A0*p0(0,t3)*A0;

path55 = p0_eq*A0*p2(0,t1)*A2*p2(2,t2)*A2*p2(2,t3)*A2;
path56 = p0_eq*A0*p2(0,t1)*A2*p2(2,t2)*A2*p1(2,t3)*A1;
path57 = p0_eq*A0*p2(0,t1)*A2*p2(2,t2)*A2*p0(2,t3)*A0;
path58 = p0_eq*A0*p2(0,t1)*A2*p1(2,t2)*A1*p2(1,t3)*A2;
path59 = p0_eq*A0*p2(0,t1)*A2*p1(2,t2)*A1*p1(1,t3)*A1;
path60 = p0_eq*A0*p2(0,t1)*A2*p1(2,t2)*A1*p0(1,t3)*A0;
path61 = p0_eq*A0*p2(0,t1)*A2*p0(2,t2)*A0*p2(0,t3)*A2;
path62 = p0_eq*A0*p2(0,t1)*A2*p0(2,t2)*A0*p1(0,t3)*A1;
path63 = p0_eq*A0*p2(0,t1)*A2*p0(2,t2)*A0*p0(0,t3)*A0;

path64 = p0_eq*A0*p1(0,t1)*A1*p2(1,t2)*A2*p2(2,t3)*A2;
path65 = p0_eq*A0*p1(0,t1)*A1*p2(1,t2)*A2*p1(2,t3)*A1;
path66 = p0_eq*A0*p1(0,t1)*A1*p2(1,t2)*A2*p0(2,t3)*A0;
path67 = p0_eq*A0*p1(0,t1)*A1*p1(1,t2)*A1*p2(1,t3)*A2;
path68 = p0_eq*A0*p1(0,t1)*A1*p1(1,t2)*A1*p1(1,t3)*A1;
path69 = p0_eq*A0*p1(0,t1)*A1*p1(1,t2)*A1*p0(1,t3)*A0;
path70 = p0_eq*A0*p1(0,t1)*A1*p0(1,t2)*A0*p2(0,t3)*A2;
path71 = p0_eq*A0*p1(0,t1)*A1*p0(1,t2)*A0*p1(0,t3)*A1;
path72 = p0_eq*A0*p1(0,t1)*A1*p0(1,t2)*A0*p0(0,t3)*A0;

path73 = p0_eq*A0*p0(0,t1)*A0*p2(0,t2)*A2*p2(2,t3)*A2;
path74 = p0_eq*A0*p0(0,t1)*A0*p2(0,t2)*A2*p1(2,t3)*A1;
path75 = p0_eq*A0*p0(0,t1)*A0*p2(0,t2)*A2*p0(2,t3)*A0;
path76 = p0_eq*A0*p0(0,t1)*A0*p1(0,t2)*A1*p2(1,t3)*A2;
path77 = p0_eq*A0*p0(0,t1)*A0*p1(0,t2)*A1*p1(1,t3)*A1;
path78 = p0_eq*A0*p0(0,t1)*A0*p1(0,t2)*A1*p0(1,t3)*A0;
path79 = p0_eq*A0*p0(0,t1)*A0*p0(0,t2)*A0*p2(0,t3)*A2;
path80 = p0_eq*A0*p0(0,t1)*A0*p0(0,t2)*A0*p1(0,t3)*A1;
path81 = p0_eq*A0*p0(0,t1)*A0*p0(0,t2)*A0*p0(0,t3)*A0;

kappa2 = path1 + path2 + path3 + path4 + path5 + path6 + path7 + path8 + path9 + path10 + path11 + path12 + path13 + path14 + path15 + path16 + path17 + path18 + path19 + path20;
kappa2 = kappa2 + path21 + path22 + path23 + path24 + path25 + path26 + path27 + path28 + path29 + path30 + path31 + path32 + path33 + path34 + path35 + path36 + path37 + path38 + path39 + path40;
kappa2 = kappa2 + path41 + path42 + path43 + path44 + path45 + path46 + path47 + path48 + path49 + path50 + path51 + path52 + path53 + path54 + path55 + path56 + path57 + path58 + path59 + path60;
kappa = kappa2 + path61 + path62 + path63 + path64 + path65 + path66 + path67 + path68 + path69 + path70 + path71 + path72 + path73 + path74 + path75 + path76 + path77 + path78 + path79 + path80 + path81;
theta_1 = (p2_eq*A2*(A2*p2(2,t1) + A1*p1(2,t1) + A0*p0(2,t1)) + p1_eq*A1*(A2*p2(1,t1) + A1*p1(1,t1) + A0*p0(1,t1)) + p0_eq*A0*(A2*p2(0,t1) + A1*p1(0,t1) + A0*p0(0,t1)));
theta_2 = (p2_eq*A2*(A2*p2(2,t3) + A1*p1(2,t3) + A0*p0(2,t3)) + p1_eq*A1*(A2*p2(1,t3) + A1*p1(1,t3) + A0*p0(1,t3)) + p0_eq*A0*(A2*p2(0,t3) + A1*p1(0,t3) + A0*p0(0,t3)));
theta = theta_1*theta_2;

tcf = kappa - theta;

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
        
        p0 = p0_eq + m*a*exp(-lam1*time) + n*x*exp(-lam2*time);
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
        
        p2 = p2_eq + m*c*exp(-lam1*time) + n*z*exp(-lam2*time);
    end

end