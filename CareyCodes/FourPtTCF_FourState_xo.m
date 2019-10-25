function [p0_eq,p1_eq,p2_eq,p3_eq,eig1,eig2,eig3] = FourPtTCF_FourState_xo(t1range,t2,A0,A1,A2,A3,k01,k10,k20,k12,k21,k23,k32)

% Detailed balance conditions
k02 = k01*k12*k20/(k10*k21);

% Define useful constants
c = (k20*k01 + k21*k01 + k21*k02)/(k20*k10 + k12*k20 + k10*k21);

% Define equilibrium constants
p0_eq = 1/(1 + c + (1+k23/k32)*(k02+k12*c)/(k20+k21));
p1_eq = c*p0_eq;
p2_eq = (k02*p0_eq + k12*p1_eq)/(k20 + k21);
p3_eq = k23/k32*p2_eq;

% Subtract the mean
Amean = p0_eq*A0 + p1_eq*A1 + p2_eq*A2 + p3_eq*A3;
A0 = A0 - Amean;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

M = [-k01-k02,k10,k20,0;k01,-k10-k12,k21,0;k02,k12,-k20-k21-k23,k32;0,0,k23,-k32];
[V,I] = eig(M);

lam = [-I(1,1),-I(2,2),-I(3,3),-I(4,4)];
[lam2,index] = sort(lam);
eig1 = lam2(2);
eig2 = lam2(3);
eig3 = lam2(4);

V2 = [V(:,index(1)),V(:,index(2)),V(:,index(3)),V(:,index(4))];

A = V2(1,2);
Ap = V2(2,2);
App = V2(3,2);
Appp = V2(4,2);
B = V2(1,3);
Bp = V2(2,3);
Bpp = V2(3,3);
Bppp = V2(4,3);
D = V2(1,4);
Dp = V2(2,4);
Dpp = V2(3,4);
Dppp = V2(4,4);

beta = (Ap*D - A*Dp)/(A*Bp - Ap*B);

% pij_t1 is prob to go from i to j in time t1
% pij(t) = pj_eq + c_1_i*v_1_j*exp(-lam1(t)) + c_2_i*v_2_j*exp(-lam2(t)*t)+ c_3_i*v_3_j*exp(-lam2(t)*t)
p00_t1 = p0_eq + c1(0)*A*exp(-lam2(2)*t1range) + c2(0)*B*exp(-lam2(3)*t1range) + c3(0)*D*exp(-lam2(4)*t1range);
p01_t1 = p1_eq + c1(0)*Ap*exp(-lam2(2)*t1range) + c2(0)*Bp*exp(-lam2(3)*t1range) + c3(0)*Dp*exp(-lam2(4)*t1range);
p02_t1 = p2_eq + c1(0)*App*exp(-lam2(2)*t1range) + c2(0)*Bpp*exp(-lam2(3)*t1range) + c3(0)*Dpp*exp(-lam2(4)*t1range);
p03_t1 = p3_eq + c1(0)*Appp*exp(-lam2(2)*t1range) + c2(0)*Bppp*exp(-lam2(3)*t1range) + c3(0)*Dppp*exp(-lam2(4)*t1range);

p00_t2 = p0_eq + c1(0)*A*exp(-lam2(2)*t2) + c2(0)*B*exp(-lam2(3)*t2) + c3(0)*D*exp(-lam2(4)*t2);
p01_t2 = p1_eq + c1(0)*Ap*exp(-lam2(2)*t2) + c2(0)*Bp*exp(-lam2(3)*t2) + c3(0)*Dp*exp(-lam2(4)*t2);
p02_t2 = p2_eq + c1(0)*App*exp(-lam2(2)*t2) + c2(0)*Bpp*exp(-lam2(3)*t2) + c3(0)*Dpp*exp(-lam2(4)*t2);
p03_t2 = p3_eq + c1(0)*Appp*exp(-lam2(2)*t2) + c2(0)*Bppp*exp(-lam2(3)*t2) + c3(0)*Dppp*exp(-lam2(4)*t2);

p00_t3 = p00_t1';
p01_t3 = p01_t1';
p02_t3 = p02_t1';
p03_t3 = p03_t1';

p10_t1 = p0_eq + c1(1)*A*exp(-lam2(2)*t1range) + c2(1)*B*exp(-lam2(3)*t1range) + c3(1)*D*exp(-lam2(4)*t1range);
p11_t1 = p1_eq + c1(1)*Ap*exp(-lam2(2)*t1range) + c2(1)*Bp*exp(-lam2(3)*t1range) + c3(1)*Dp*exp(-lam2(4)*t1range);
p12_t1 = p2_eq + c1(1)*App*exp(-lam2(2)*t1range) + c2(1)*Bpp*exp(-lam2(3)*t1range) + c3(1)*Dpp*exp(-lam2(4)*t1range);
p13_t1 = p3_eq + c1(1)*Appp*exp(-lam2(2)*t1range) + c2(1)*Bppp*exp(-lam2(3)*t1range) + c3(1)*Dppp*exp(-lam2(4)*t1range);

p10_t2 = p0_eq + c1(1)*A*exp(-lam2(2)*t2) + c2(1)*B*exp(-lam2(3)*t2) + c3(1)*D*exp(-lam2(4)*t2);
p11_t2 = p1_eq + c1(1)*Ap*exp(-lam2(2)*t2) + c2(1)*Bp*exp(-lam2(3)*t2) + c3(1)*Dp*exp(-lam2(4)*t2);
p12_t2 = p2_eq + c1(1)*App*exp(-lam2(2)*t2) + c2(1)*Bpp*exp(-lam2(3)*t2) + c3(1)*Dpp*exp(-lam2(4)*t2);
p13_t2 = p3_eq + c1(1)*Appp*exp(-lam2(2)*t2) + c2(1)*Bppp*exp(-lam2(3)*t2) + c3(1)*Dppp*exp(-lam2(4)*t2);

p10_t3 = p10_t1';
p11_t3 = p11_t1';
p12_t3 = p12_t1';
p13_t3 = p13_t1';

p20_t1 = p0_eq + c1(2)*A*exp(-lam2(2)*t1range) + c2(2)*B*exp(-lam2(3)*t1range) + c3(2)*D*exp(-lam2(4)*t1range);
p21_t1 = p1_eq + c1(2)*Ap*exp(-lam2(2)*t1range) + c2(2)*Bp*exp(-lam2(3)*t1range) + c3(2)*Dp*exp(-lam2(4)*t1range);
p22_t1 = p2_eq + c1(2)*App*exp(-lam2(2)*t1range) + c2(2)*Bpp*exp(-lam2(3)*t1range) + c3(2)*Dpp*exp(-lam2(4)*t1range);
p23_t1 = p3_eq + c1(2)*Appp*exp(-lam2(2)*t1range) + c2(2)*Bppp*exp(-lam2(3)*t1range) + c3(2)*Dppp*exp(-lam2(4)*t1range);

p20_t2 = p0_eq + c1(2)*A*exp(-lam2(2)*t2) + c2(2)*B*exp(-lam2(3)*t2) + c3(2)*D*exp(-lam2(4)*t2);
p21_t2 = p1_eq + c1(2)*Ap*exp(-lam2(2)*t2) + c2(2)*Bp*exp(-lam2(3)*t2) + c3(2)*Dp*exp(-lam2(4)*t2);
p22_t2 = p2_eq + c1(2)*App*exp(-lam2(2)*t2) + c2(2)*Bpp*exp(-lam2(3)*t2) + c3(2)*Dpp*exp(-lam2(4)*t2);
p23_t2 = p3_eq + c1(2)*Appp*exp(-lam2(2)*t2) + c2(2)*Bppp*exp(-lam2(3)*t2) + c3(2)*Dppp*exp(-lam2(4)*t2);

p20_t3 = p20_t1';
p21_t3 = p21_t1';
p22_t3 = p22_t1';
p23_t3 = p23_t1';

p30_t1 = p0_eq + c1(3)*A*exp(-lam2(2)*t1range) + c2(3)*B*exp(-lam2(3)*t1range) + c3(3)*D*exp(-lam2(4)*t1range);
p31_t1 = p1_eq + c1(3)*Ap*exp(-lam2(2)*t1range) + c2(3)*Bp*exp(-lam2(3)*t1range) + c3(3)*Dp*exp(-lam2(4)*t1range);
p32_t1 = p2_eq + c1(3)*App*exp(-lam2(2)*t1range) + c2(3)*Bpp*exp(-lam2(3)*t1range) + c3(3)*Dpp*exp(-lam2(4)*t1range);
p33_t1 = p3_eq + c1(3)*Appp*exp(-lam2(2)*t1range) + c2(3)*Bppp*exp(-lam2(3)*t1range) + c3(3)*Dppp*exp(-lam2(4)*t1range);

p30_t2 = p0_eq + c1(3)*A*exp(-lam2(2)*t2) + c2(3)*B*exp(-lam2(3)*t2) + c3(3)*D*exp(-lam2(4)*t2);
p31_t2 = p1_eq + c1(3)*Ap*exp(-lam2(2)*t2) + c2(3)*Bp*exp(-lam2(3)*t2) + c3(3)*Dp*exp(-lam2(4)*t2);
p32_t2 = p2_eq + c1(3)*App*exp(-lam2(2)*t2) + c2(3)*Bpp*exp(-lam2(3)*t2) + c3(3)*Dpp*exp(-lam2(4)*t2);
p33_t2 = p3_eq + c1(3)*Appp*exp(-lam2(2)*t2) + c2(3)*Bppp*exp(-lam2(3)*t2) + c3(3)*Dppp*exp(-lam2(4)*t2);

p30_t3 = p30_t1';
p31_t3 = p31_t1';
p32_t3 = p32_t1';
p33_t3 = p33_t1';

% Define all paths
kappa1 = path(0,0,0,0) + path(0,0,0,1) + path(0,0,0,2) + path(0,0,0,3) + path(0,0,1,0) + path(0,0,1,1) + path(0,0,1,2) + path(0,0,1,3);
kappa2 = path(0,0,2,0) + path(0,0,2,1) + path(0,0,2,2) + path(0,0,2,3) + path(0,0,3,0) + path(0,0,3,1) + path(0,0,3,2) + path(0,0,3,3);
kappa3 = path(0,1,0,0) + path(0,1,0,1) + path(0,1,0,2) + path(0,1,0,3) + path(0,1,1,0) + path(0,1,1,1) + path(0,1,1,2) + path(0,1,1,3);
kappa4 = path(0,1,2,0) + path(0,1,2,1) + path(0,1,2,2) + path(0,1,2,3) + path(0,1,3,0) + path(0,1,3,1) + path(0,1,3,2) + path(0,1,3,3);
kappa5 = path(0,2,0,0) + path(0,2,0,1) + path(0,2,0,2) + path(0,2,0,3) + path(0,2,1,0) + path(0,2,1,1) + path(0,2,1,2) + path(0,2,1,3);
kappa6 = path(0,2,2,0) + path(0,2,2,1) + path(0,2,2,2) + path(0,2,2,3) + path(0,2,3,0) + path(0,2,3,1) + path(0,2,3,2) + path(0,2,3,3);
kappa7 = path(0,3,0,0) + path(0,3,0,1) + path(0,3,0,2) + path(0,3,0,3) + path(0,3,1,0) + path(0,3,1,1) + path(0,3,1,2) + path(0,3,1,3);
kappa8 = path(0,3,2,0) + path(0,3,2,1) + path(0,3,2,2) + path(0,3,2,3) + path(0,3,3,0) + path(0,3,3,1) + path(0,3,3,2) + path(0,3,3,3);

kappa9 = path(1,0,0,0) + path(1,0,0,1) + path(1,0,0,2) + path(1,0,0,3) + path(1,0,1,0) + path(1,0,1,1) + path(1,0,1,2) + path(1,0,1,3);
kappa10 = path(1,0,2,0) + path(1,0,2,1) + path(1,0,2,2) + path(1,0,2,3) + path(1,0,3,0) + path(1,0,3,1) + path(1,0,3,2) + path(1,0,3,3);
kappa11 = path(1,1,0,0) + path(1,1,0,1) + path(1,1,0,2) + path(1,1,0,3) + path(1,1,1,0) + path(1,1,1,1) + path(1,1,1,2) + path(1,1,1,3);
kappa12 = path(1,1,2,0) + path(1,1,2,1) + path(1,1,2,2) + path(1,1,2,3) + path(1,1,3,0) + path(1,1,3,1) + path(1,1,3,2) + path(1,1,3,3);
kappa13 = path(1,2,0,0) + path(1,2,0,1) + path(1,2,0,2) + path(1,2,0,3) + path(1,2,1,0) + path(1,2,1,1) + path(1,2,1,2) + path(1,2,1,3);
kappa14 = path(1,2,2,0) + path(1,2,2,1) + path(1,2,2,2) + path(1,2,2,3) + path(1,2,3,0) + path(1,2,3,1) + path(1,2,3,2) + path(1,2,3,3);
kappa15 = path(1,3,0,0) + path(1,3,0,1) + path(1,3,0,2) + path(1,3,0,3) + path(1,3,1,0) + path(1,3,1,1) + path(1,3,1,2) + path(1,3,1,3);
kappa16 = path(1,3,2,0) + path(1,3,2,1) + path(1,3,2,2) + path(1,3,2,3) + path(1,3,3,0) + path(1,3,3,1) + path(1,3,3,2) + path(1,3,3,3);

kappa17 = path(2,0,0,0) + path(2,0,0,1) + path(2,0,0,2) + path(2,0,0,3) + path(2,0,1,0) + path(2,0,1,1) + path(2,0,1,2) + path(2,0,1,3);
kappa18 = path(2,0,2,0) + path(2,0,2,1) + path(2,0,2,2) + path(2,0,2,3) + path(2,0,3,0) + path(2,0,3,1) + path(2,0,3,2) + path(2,0,3,3);
kappa19 = path(2,1,0,0) + path(2,1,0,1) + path(2,1,0,2) + path(2,1,0,3) + path(2,1,1,0) + path(2,1,1,1) + path(2,1,1,2) + path(2,1,1,3);
kappa20 = path(2,1,2,0) + path(2,1,2,1) + path(2,1,2,2) + path(2,1,2,3) + path(2,1,3,0) + path(2,1,3,1) + path(2,1,3,2) + path(2,1,3,3);
kappa21 = path(2,2,0,0) + path(2,2,0,1) + path(2,2,0,2) + path(2,2,0,3) + path(2,2,1,0) + path(2,2,1,1) + path(2,2,1,2) + path(2,2,1,3);
kappa22 = path(2,2,2,0) + path(2,2,2,1) + path(2,2,2,2) + path(2,2,2,3) + path(2,2,3,0) + path(2,2,3,1) + path(2,2,3,2) + path(2,2,3,3);
kappa23 = path(2,3,0,0) + path(2,3,0,1) + path(2,3,0,2) + path(2,3,0,3) + path(2,3,1,0) + path(2,3,1,1) + path(2,3,1,2) + path(2,3,1,3);
kappa24 = path(2,3,2,0) + path(2,3,2,1) + path(2,3,2,2) + path(2,3,2,3) + path(2,3,3,0) + path(2,3,3,1) + path(2,3,3,2) + path(2,3,3,3);

kappa25 = path(3,0,0,0) + path(3,0,0,1) + path(3,0,0,2) + path(3,0,0,3) + path(3,0,1,0) + path(3,0,1,1) + path(3,0,1,2) + path(3,0,1,3);
kappa26 = path(3,0,2,0) + path(3,0,2,1) + path(3,0,2,2) + path(3,0,2,3) + path(3,0,3,0) + path(3,0,3,1) + path(3,0,3,2) + path(3,0,3,3);
kappa27 = path(3,1,0,0) + path(3,1,0,1) + path(3,1,0,2) + path(3,1,0,3) + path(3,1,1,0) + path(3,1,1,1) + path(3,1,1,2) + path(3,1,1,3);
kappa28 = path(3,1,2,0) + path(3,1,2,1) + path(3,1,2,2) + path(3,1,2,3) + path(3,1,3,0) + path(3,1,3,1) + path(3,1,3,2) + path(3,1,3,3);
kappa29 = path(3,2,0,0) + path(3,2,0,1) + path(3,2,0,2) + path(3,2,0,3) + path(3,2,1,0) + path(3,2,1,1) + path(3,2,1,2) + path(3,2,1,3);
kappa30 = path(3,2,2,0) + path(3,2,2,1) + path(3,2,2,2) + path(3,2,2,3) + path(3,2,3,0) + path(3,2,3,1) + path(3,2,3,2) + path(3,2,3,3);
kappa31 = path(3,3,0,0) + path(3,3,0,1) + path(3,3,0,2) + path(3,3,0,3) + path(3,3,1,0) + path(3,3,1,1) + path(3,3,1,2) + path(3,3,1,3);
kappa32 = path(3,3,2,0) + path(3,3,2,1) + path(3,3,2,2) + path(3,3,2,3) + path(3,3,3,0) + path(3,3,3,1) + path(3,3,3,2) + path(3,3,3,3);

kappa = kappa1 + kappa2 + kappa3 + kappa4 + kappa5 + kappa6 + kappa7 + kappa8 + kappa9 + kappa10;
kappa = kappa + kappa11 + kappa12 + kappa13 + kappa14 + kappa15 + kappa16 + kappa17 + kappa18 + kappa19 + kappa20;
kappa = kappa + kappa21 + kappa22 + kappa23 + kappa24 + kappa25 + kappa26 + kappa27 + kappa28 + kappa29 + kappa30 + kappa31 + kappa32;

theta1 = path2(0,0,0,0) + path2(0,0,0,1) + path2(0,0,0,2) + path2(0,0,0,3) + path2(0,0,1,0) + path2(0,0,1,1) + path2(0,0,1,2) + path2(0,0,1,3);
theta2 = path2(0,0,2,0) + path2(0,0,2,1) + path2(0,0,2,2) + path2(0,0,2,3) + path2(0,0,3,0) + path2(0,0,3,1) + path2(0,0,3,2) + path2(0,0,3,3);
theta3 = path2(0,1,0,0) + path2(0,1,0,1) + path2(0,1,0,2) + path2(0,1,0,3) + path2(0,1,1,0) + path2(0,1,1,1) + path2(0,1,1,2) + path2(0,1,1,3);
theta4 = path2(0,1,2,0) + path2(0,1,2,1) + path2(0,1,2,2) + path2(0,1,2,3) + path2(0,1,3,0) + path2(0,1,3,1) + path2(0,1,3,2) + path2(0,1,3,3);
theta5 = path2(0,2,0,0) + path2(0,2,0,1) + path2(0,2,0,2) + path2(0,2,0,3) + path2(0,2,1,0) + path2(0,2,1,1) + path2(0,2,1,2) + path2(0,2,1,3);
theta6 = path2(0,2,2,0) + path2(0,2,2,1) + path2(0,2,2,2) + path2(0,2,2,3) + path2(0,2,3,0) + path2(0,2,3,1) + path2(0,2,3,2) + path2(0,2,3,3);
theta7 = path2(0,3,0,0) + path2(0,3,0,1) + path2(0,3,0,2) + path2(0,3,0,3) + path2(0,3,1,0) + path2(0,3,1,1) + path2(0,3,1,2) + path2(0,3,1,3);
theta8 = path2(0,3,2,0) + path2(0,3,2,1) + path2(0,3,2,2) + path2(0,3,2,3) + path2(0,3,3,0) + path2(0,3,3,1) + path2(0,3,3,2) + path2(0,3,3,3);

theta9 = path2(1,0,0,0) + path2(1,0,0,1) + path2(1,0,0,2) + path2(1,0,0,3) + path2(1,0,1,0) + path2(1,0,1,1) + path2(1,0,1,2) + path2(1,0,1,3);
theta10 = path2(1,0,2,0) + path2(1,0,2,1) + path2(1,0,2,2) + path2(1,0,2,3) + path2(1,0,3,0) + path2(1,0,3,1) + path2(1,0,3,2) + path2(1,0,3,3);
theta11 = path2(1,1,0,0) + path2(1,1,0,1) + path2(1,1,0,2) + path2(1,1,0,3) + path2(1,1,1,0) + path2(1,1,1,1) + path2(1,1,1,2) + path2(1,1,1,3);
theta12 = path2(1,1,2,0) + path2(1,1,2,1) + path2(1,1,2,2) + path2(1,1,2,3) + path2(1,1,3,0) + path2(1,1,3,1) + path2(1,1,3,2) + path2(1,1,3,3);
theta13 = path2(1,2,0,0) + path2(1,2,0,1) + path2(1,2,0,2) + path2(1,2,0,3) + path2(1,2,1,0) + path2(1,2,1,1) + path2(1,2,1,2) + path2(1,2,1,3);
theta14 = path2(1,2,2,0) + path2(1,2,2,1) + path2(1,2,2,2) + path2(1,2,2,3) + path2(1,2,3,0) + path2(1,2,3,1) + path2(1,2,3,2) + path2(1,2,3,3);
theta15 = path2(1,3,0,0) + path2(1,3,0,1) + path2(1,3,0,2) + path2(1,3,0,3) + path2(1,3,1,0) + path2(1,3,1,1) + path2(1,3,1,2) + path2(1,3,1,3);
theta16 = path2(1,3,2,0) + path2(1,3,2,1) + path2(1,3,2,2) + path2(1,3,2,3) + path2(1,3,3,0) + path2(1,3,3,1) + path2(1,3,3,2) + path2(1,3,3,3);

theta17 = path2(2,0,0,0) + path2(2,0,0,1) + path2(2,0,0,2) + path2(2,0,0,3) + path2(2,0,1,0) + path2(2,0,1,1) + path2(2,0,1,2) + path2(2,0,1,3);
theta18 = path2(2,0,2,0) + path2(2,0,2,1) + path2(2,0,2,2) + path2(2,0,2,3) + path2(2,0,3,0) + path2(2,0,3,1) + path2(2,0,3,2) + path2(2,0,3,3);
theta19 = path2(2,1,0,0) + path2(2,1,0,1) + path2(2,1,0,2) + path2(2,1,0,3) + path2(2,1,1,0) + path2(2,1,1,1) + path2(2,1,1,2) + path2(2,1,1,3);
theta20 = path2(2,1,2,0) + path2(2,1,2,1) + path2(2,1,2,2) + path2(2,1,2,3) + path2(2,1,3,0) + path2(2,1,3,1) + path2(2,1,3,2) + path2(2,1,3,3);
theta21 = path2(2,2,0,0) + path2(2,2,0,1) + path2(2,2,0,2) + path2(2,2,0,3) + path2(2,2,1,0) + path2(2,2,1,1) + path2(2,2,1,2) + path2(2,2,1,3);
theta22 = path2(2,2,2,0) + path2(2,2,2,1) + path2(2,2,2,2) + path2(2,2,2,3) + path2(2,2,3,0) + path2(2,2,3,1) + path2(2,2,3,2) + path2(2,2,3,3);
theta23 = path2(2,3,0,0) + path2(2,3,0,1) + path2(2,3,0,2) + path2(2,3,0,3) + path2(2,3,1,0) + path2(2,3,1,1) + path2(2,3,1,2) + path2(2,3,1,3);
theta24 = path2(2,3,2,0) + path2(2,3,2,1) + path2(2,3,2,2) + path2(2,3,2,3) + path2(2,3,3,0) + path2(2,3,3,1) + path2(2,3,3,2) + path2(2,3,3,3);

theta25 = path2(3,0,0,0) + path2(3,0,0,1) + path2(3,0,0,2) + path2(3,0,0,3) + path2(3,0,1,0) + path2(3,0,1,1) + path2(3,0,1,2) + path2(3,0,1,3);
theta26 = path2(3,0,2,0) + path2(3,0,2,1) + path2(3,0,2,2) + path2(3,0,2,3) + path2(3,0,3,0) + path2(3,0,3,1) + path2(3,0,3,2) + path2(3,0,3,3);
theta27 = path2(3,1,0,0) + path2(3,1,0,1) + path2(3,1,0,2) + path2(3,1,0,3) + path2(3,1,1,0) + path2(3,1,1,1) + path2(3,1,1,2) + path2(3,1,1,3);
theta28 = path2(3,1,2,0) + path2(3,1,2,1) + path2(3,1,2,2) + path2(3,1,2,3) + path2(3,1,3,0) + path2(3,1,3,1) + path2(3,1,3,2) + path2(3,1,3,3);
theta29 = path2(3,2,0,0) + path2(3,2,0,1) + path2(3,2,0,2) + path2(3,2,0,3) + path2(3,2,1,0) + path2(3,2,1,1) + path2(3,2,1,2) + path2(3,2,1,3);
theta30 = path2(3,2,2,0) + path2(3,2,2,1) + path2(3,2,2,2) + path2(3,2,2,3) + path2(3,2,3,0) + path2(3,2,3,1) + path2(3,2,3,2) + path2(3,2,3,3);
theta31 = path2(3,3,0,0) + path2(3,3,0,1) + path2(3,3,0,2) + path2(3,3,0,3) + path2(3,3,1,0) + path2(3,3,1,1) + path2(3,3,1,2) + path2(3,3,1,3);
theta32 = path2(3,3,2,0) + path2(3,3,2,1) + path2(3,3,2,2) + path2(3,3,2,3) + path2(3,3,3,0) + path2(3,3,3,1) + path2(3,3,3,2) + path2(3,3,3,3);

theta = theta1 + theta2 + theta3 + theta4 + theta5 + theta6 + theta7 + theta8 + theta9 + theta10;
theta = theta + theta11 + theta12 + theta13 + theta14 + theta15 + theta16 + theta17 + theta18 + theta19 + theta20;
theta = theta + theta21 + theta22 + theta23 + theta24 + theta25 + theta26 + theta27 + theta28 + theta29 + theta30 + theta31 + theta32;

tcf = kappa - theta;

    function out = alpha(st)
            if st == 0
                out = (A*(-p1_eq) - Ap*(1 - p0_eq))/(A*Bp - Ap*B);
            elseif st == 1
                out = (A*(1-p1_eq) - Ap*(-p0_eq))/(A*Bp - Ap*B);
            else
                out = (A*(-p1_eq) - Ap*(-p0_eq))/(A*Bp - Ap*B);
            end
    end
    
    function out = c3(st)
        if st == 0
            out = (A*(-p2_eq - Bpp*alpha(0)) - App*(1-p0_eq - B*alpha(0)))/(A*Bpp*beta + A*Dpp - beta*B*App - D*App);
        elseif st == 1
            out = (A*(-p2_eq - Bpp*alpha(1)) - App*(-p0_eq - B*alpha(1)))/(A*Bpp*beta + A*Dpp - beta*B*App - D*App);
        elseif st == 2
            out = (A*(1 - p2_eq - Bpp*alpha(2)) - App*(-p0_eq - B*alpha(2)))/(A*Bpp*beta + A*Dpp - beta*B*App - D*App);
        else
            out = (A*(-p2_eq - Bpp*alpha(3)) - App*(-p0_eq - B*alpha(3)))/(A*Bpp*beta + A*Dpp - beta*B*App - D*App);
        end
    end

    function out = c2(st)
        if st == 0
            out = alpha(0) + beta*c3(0);
        elseif st == 1
            out = alpha(1) + beta*c3(1);
        elseif st == 2
            out = alpha(2) + beta*c3(2);
        else
            out = alpha(3) + beta*c3(3);
        end
    end

    function out = c1(st)
        if st == 0
            out = (1 - p0_eq - c2(0)*B - c3(0)*D)/A;
        elseif st == 1
            out = (- p0_eq - c2(1)*B - c3(1)*D)/A;
        elseif st == 2
            out = (- p0_eq - c2(2)*B - c3(2)*D)/A;
        else
            out = (- p0_eq - c2(3)*B - c3(3)*D)/A;
        end
    end

function out = path(st1,st2,st3,st4)
    if st1 == 0
        peq = p0_eq;
        Ai = A0;
        if st2 == 0
            pij = p00_t1;
        elseif st2 == 1
            pij = p01_t1;
        elseif st2 == 2
            pij = p02_t1;
        else
            pij = p03_t1;
        end
    elseif st1 == 1
        peq = p1_eq;
        Ai = A1;
        if st2 == 0
            pij = p10_t1;
        elseif st2 == 1
            pij = p11_t1;
        elseif st2 == 2
            pij = p12_t1;
        else
            pij = p13_t1;
        end
    elseif st1 == 2
        peq = p2_eq;
        Ai = A2;
        if st2 == 0
            pij = p20_t1;
        elseif st2 == 1
            pij = p21_t1;
        elseif st2 == 2
            pij = p22_t1;
        else
            pij = p23_t1;
        end
    else
        if st2 == 0
            pij = p30_t1;
        elseif st2 == 1
            pij = p31_t1;
        elseif st2 == 2
            pij = p32_t1;
        else
            pij = p33_t1;
        end
        peq = p3_eq;
        Ai = A3;
    end
    
    if st2 == 0
        Aj = A0;
        if st3 == 0
            pjk = p00_t2;
        elseif st3 == 1
            pjk = p01_t2;
        elseif st3 == 2
            pjk = p02_t2;
        else
            pjk = p03_t2;
        end
    elseif st2 == 1
        Aj = A1;
        if st3 == 0
            pjk = p10_t2;
        elseif st3 == 1
            pjk = p11_t2;
        elseif st3 == 2
            pjk = p12_t2;
        else
            pjk = p13_t2;
        end
    elseif st2 == 2
        Aj = A2;
        if st3 == 0
            pjk = p20_t2;
        elseif st3 == 1
            pjk = p21_t2;
        elseif st3 == 2
            pjk = p22_t2;
        else
            pjk = p23_t2;
        end
    else
        Aj = A3;
        if st3 == 0
            pjk = p30_t2;
        elseif st3 == 1
            pjk = p31_t2;
        elseif st3 == 2
            pjk = p32_t2;
        else
            pjk = p33_t2;
        end
    end
    
    if st3 == 0
        Ak = A0;
        peq2 = p0_eq;
        if st4 == 0
            pkl = p00_t3;
        elseif st4 == 1
            pkl = p01_t3;
        elseif st4 == 2
            pkl = p02_t3;
        else
            pkl = p03_t3;
        end
    elseif st3 == 1
        Ak = A1;
        peq2 = p1_eq;
        if st4 == 0
            pkl = p10_t3;
        elseif st4 == 1
            pkl = p11_t3;
        elseif st4 == 2
            pkl = p12_t3;
        else
            pkl = p13_t3;
        end
    elseif st3 == 2
        Ak = A2;
        peq2 = p2_eq;
        if st4 == 0
            pkl = p20_t3;
        elseif st4 == 1
            pkl = p21_t3;
        elseif st4 == 2
            pkl = p22_t3;
        else
            pkl = p23_t3;
        end
    else
        Ak = A3;
        peq2 = p3_eq;
        if st4 == 0
            pkl = p30_t3;
        elseif st4 == 1
            pkl = p31_t3;
        elseif st4 == 2
            pkl = p32_t3;
        else
            pkl = p33_t3;
        end
    end
    
    if st4 == 0
        Al = A0;
    elseif st4 == 1
        Al = A1;
    elseif st4 == 2
        Al = A2;
    else
        Al = A3;
    end
    
    out = Ai*Aj*Ak*Al*peq*pij*pjk*pkl;
    
end

function out = path2(st1,st2,st3,st4)
    if st1 == 0
        peq = p0_eq;
        Ai = A0;
        if st2 == 0
            pij = p00_t1;
        elseif st2 == 1
            pij = p01_t1;
        elseif st2 == 2
            pij = p02_t1;
        else
            pij = p03_t1;
        end
    elseif st1 == 1
        peq = p1_eq;
        Ai = A1;
        if st2 == 0
            pij = p10_t1;
        elseif st2 == 1
            pij = p11_t1;
        elseif st2 == 2
            pij = p12_t1;
        else
            pij = p13_t1;
        end
    elseif st1 == 2
        peq = p2_eq;
        Ai = A2;
        if st2 == 0
            pij = p20_t1;
        elseif st2 == 1
            pij = p21_t1;
        elseif st2 == 2
            pij = p22_t1;
        else
            pij = p23_t1;
        end
    else
        if st2 == 0
            pij = p30_t1;
        elseif st2 == 1
            pij = p31_t1;
        elseif st2 == 2
            pij = p32_t1;
        else
            pij = p33_t1;
        end
        peq = p3_eq;
        Ai = A3;
    end
    
    if st2 == 0
        Aj = A0;
        if st3 == 0
            pjk = p00_t2;
        elseif st3 == 1
            pjk = p01_t2;
        elseif st3 == 2
            pjk = p02_t2;
        else
            pjk = p03_t2;
        end
    elseif st2 == 1
        Aj = A1;
        if st3 == 0
            pjk = p10_t2;
        elseif st3 == 1
            pjk = p11_t2;
        elseif st3 == 2
            pjk = p12_t2;
        else
            pjk = p13_t2;
        end
    elseif st2 == 2
        Aj = A2;
        if st3 == 0
            pjk = p20_t2;
        elseif st3 == 1
            pjk = p21_t2;
        elseif st3 == 2
            pjk = p22_t2;
        else
            pjk = p23_t2;
        end
    else
        Aj = A3;
        if st3 == 0
            pjk = p30_t2;
        elseif st3 == 1
            pjk = p31_t2;
        elseif st3 == 2
            pjk = p32_t2;
        else
            pjk = p33_t2;
        end
    end
    
    if st3 == 0
        Ak = A0;
        peq2 = p0_eq;
        if st4 == 0
            pkl = p00_t3;
        elseif st4 == 1
            pkl = p01_t3;
        elseif st4 == 2
            pkl = p02_t3;
        else
            pkl = p03_t3;
        end
    elseif st3 == 1
        Ak = A1;
        peq2 = p1_eq;
        if st4 == 0
            pkl = p10_t3;
        elseif st4 == 1
            pkl = p11_t3;
        elseif st4 == 2
            pkl = p12_t3;
        else
            pkl = p13_t3;
        end
    elseif st3 == 2
        Ak = A2;
        peq2 = p2_eq;
        if st4 == 0
            pkl = p20_t3;
        elseif st4 == 1
            pkl = p21_t3;
        elseif st4 == 2
            pkl = p22_t3;
        else
            pkl = p23_t3;
        end
    else
        Ak = A3;
        peq2 = p3_eq;
        if st4 == 0
            pkl = p30_t3;
        elseif st4 == 1
            pkl = p31_t3;
        elseif st4 == 2
            pkl = p32_t3;
        else
            pkl = p33_t3;
        end
    end
    
    if st4 == 0
        Al = A0;
    elseif st4 == 1
        Al = A1;
    elseif st4 == 2
        Al = A2;
    else
        Al = A3;
    end
    
    out = Ai*Aj*Ak*Al*peq*pij*pkl*peq2;
    
end
end