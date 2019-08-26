function tcf = tcf_2state(time,A0,A1,k01,k10)

% Define eigenvalue lam
lam = k01 + k10;

p1_eq = 1/(k10/k01 + 1);
p0_eq = k10/k01*p1_eq;

p00 = p0_eq + (1-p0_eq)*exp(-lam*time);
p01 = p1_eq - p1_eq*exp(-lam*time);
p10 = p0_eq - p0_eq*exp(-lam*time);
p11 = p1_eq + (1-p1_eq)*exp(-lam*time);

tcf = A0*A0*p0_eq*p00 + A0*A1*p0_eq*p01 + A1*A0*p1_eq*p10 + A1*A1*p1_eq*p11;

semilogx(time,tcf)

end