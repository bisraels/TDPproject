    
syms k01 k02 k10 k20 k12 k13 k14 k15 k16 k21 k23 k24 k25 k26 k31 k32 k34 k35 k36 k41 k42 k43 k45 k46 k51 k52 k53 k54 k56 k61 k62 k63 k64 k65 
syms P1 P2 P3 P3 P4


K = [-(k01 + k02 + 0), k10, k20, 0;...
    k01, -(k10 + k12 + 0), k21, 0;...
    k02, k12, -(k20 + k21 + k23), k32;...
    0, 0, k23, -(0 + 0 + k32);];

    P = [sym('P',[1 4],'real')];
    
    eqns = [mtimes(K,P') == 0; P1 + P2 + P3 + P4 == 1];
    
    solns = solve(eqns,[P1 P2 P3 P4]);
    P1_n4 = solns.P1
    P2_n4 = solns.P2
    P3_n4 = solns.P3
    P4_n4 = solns.P4
    
    %Note: use subs to evaluate a symbolic expression with its values