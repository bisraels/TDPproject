function [kappa,theta,tcf] = FourPtTCF_2state(tau1range,tau2,tau3range,A0,A1,k01,k10)
global k10 k01 p1_eq p0_eq

tcf = zeros(length(tau1range),length(tau3range));
kappa = zeros(length(tau1range),length(tau3range));
theta = zeros(length(tau1range),length(tau3range));

p1_eq = 1/(1 + k10/k01);
p0_eq = k10/k01*p1_eq;

for k = 1:length(tau1range)
    t1 = tau1range(k);
    t2 = tau2;
    t3 = tau3range;
    
path1 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1*p1(1,t3)*A1;
path2 = p1_eq*A1*p1(1,t1)*A1*p1(1,t2)*A1*p0(0,t3)*A0;
path3 = p1_eq*A1*p1(1,t1)*A1*p0(0,t2)*A0*p1(0,t3)*A1;
path4 = p1_eq*A1*p1(1,t1)*A1*p0(0,t2)*A0*p0(1,t3)*A0;
path5 = p1_eq*A1*p0(0,t1)*A0*p1(0,t2)*A1*p1(1,t3)*A1;
path6 = p1_eq*A1*p0(0,t1)*A0*p1(0,t2)*A1*p0(0,t3)*A0;
path7 = p1_eq*A1*p0(0,t1)*A0*p0(1,t2)*A0*p1(0,t3)*A1;
path8 = p1_eq*A1*p0(0,t1)*A0*p0(1,t2)*A0*p0(1,t3)*A0;
path9 = p0_eq*A0*p1(0,t1)*A1*p1(1,t2)*A1*p1(1,t3)*A1;
path10 = p0_eq*A0*p1(0,t1)*A1*p1(1,t2)*A1*p0(0,t3)*A0;
path11 = p0_eq*A0*p1(0,t1)*A1*p0(0,t2)*A0*p1(0,t3)*A1;
path12 = p0_eq*A0*p1(0,t1)*A1*p0(0,t2)*A0*p0(1,t3)*A0;
path13 = p0_eq*A0*p0(1,t1)*A0*p1(0,t2)*A1*p1(1,t3)*A1;
path14 = p0_eq*A0*p0(1,t1)*A0*p1(0,t2)*A1*p0(0,t3)*A0;
path15 = p0_eq*A0*p0(1,t1)*A0*p0(1,t2)*A0*p1(0,t3)*A1;
path16 = p0_eq*A0*p0(1,t1)*A0*p0(1,t2)*A0*p0(1,t3)*A0;

kappa(k,:) = path1 + path2 + path3 + path4 + path5 + path6 + path7 + path8 + path9 + path10 + path11 + path12 + path13 + path14 + path15 + path16;
theta(k,:) = (p1_eq*A1*(A1*p1(1,t1) + A0*p0(0,t1)) + p0_eq*A0*(A1*p1(0,t1) + A0*p0(1,t1)))*(p1_eq*A1*(A1*p1(1,t3) + A0*p0(0,t3)) + p0_eq*A0*(A1*p1(0,t3) + A0*p0(1,t3)));

end

tcf = kappa - theta;
end

function p1 = p1(p1_init,time)
global k10 k01 p1_eq p0_eq
    p1 = p1_eq + (p1_init - p1_eq)*exp(-(k01+k10)*time);
end

function p0 = p0(p0_init,time)
global k10 k01 p1_eq p0_eq
    p0 = p0_eq + (p0_init - p0_eq)*exp(-(k01+k10)*time);
end