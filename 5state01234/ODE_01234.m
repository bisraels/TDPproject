% ODE_01234



 function dPdt = ODE_01234(t,P,opts,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43)
%function dPdt = ODE_01234(t,P,opt,K)
dPdt = zeros(5,1);
K= [-(k01+k02+k03+k04), k10,             k20,              k30,              k40;...
        k01,      -(k10+k12+k13+k14),     k21,              k31,              k41;...
        k02,             k12,       -(k20+k21+k23+k24),     k32,              k42;...
        k03,             k13,             k23,      -(k30+k31+k32+k34),       k43;...
        k04,             k14,             k24,              k34,        -(k40+k41+k42+k43)];
KP = K*P;
dPdt(1) = KP(1);
dPdt(2) = KP(2);
dPdt(3) = KP(3);
dPdt(4) = KP(4);
dPdt(5) = KP(5);

end