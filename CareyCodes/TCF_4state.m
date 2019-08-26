
% Calculates 2-pt. TCF using analytical formula
% Is called by ..... fmincon_gp32_4state_1p0uMgp32
function tcf = TCF_4state(time,A0,A1,A2,A3,k01,k10,k20,k12,k21,k23,k32)

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

p00 = p0_eq + c1(0)*A*exp(-lam2(2)*time) + c2(0)*B*exp(-lam2(3)*time) + c3(0)*D*exp(-lam2(4)*time);
p01 = p1_eq + c1(0)*Ap*exp(-lam2(2)*time) + c2(0)*Bp*exp(-lam2(3)*time) + c3(0)*Dp*exp(-lam2(4)*time);
p02 = p2_eq + c1(0)*App*exp(-lam2(2)*time) + c2(0)*Bpp*exp(-lam2(3)*time) + c3(0)*Dpp*exp(-lam2(4)*time);
p03 = p3_eq + c1(0)*Appp*exp(-lam2(2)*time) + c2(0)*Bppp*exp(-lam2(3)*time) + c3(0)*Dppp*exp(-lam2(4)*time);

p10 = p0_eq + c1(1)*A*exp(-lam2(2)*time) + c2(1)*V2(1,3)*exp(-lam2(3)*time) + c3(1)*V2(1,4)*exp(-lam2(4)*time);
p11 = p1_eq + c1(1)*Ap*exp(-lam2(2)*time) + c2(1)*V2(2,3)*exp(-lam2(3)*time) + c3(1)*V2(2,4)*exp(-lam2(4)*time);
p12 = p2_eq + c1(1)*App*exp(-lam2(2)*time) + c2(1)*V2(3,3)*exp(-lam2(3)*time) + c3(1)*V2(3,4)*exp(-lam2(4)*time);
p13 = p3_eq + c1(1)*Appp*exp(-lam2(2)*time) + c2(1)*V2(4,3)*exp(-lam2(3)*time) + c3(1)*V2(4,4)*exp(-lam2(4)*time);

p20 = p0_eq + c1(2)*A*exp(-lam2(2)*time) + c2(2)*B*exp(-lam2(3)*time) + c3(2)*D*exp(-lam2(4)*time);
p21 = p1_eq + c1(2)*Ap*exp(-lam2(2)*time) + c2(2)*Bp*exp(-lam2(3)*time) + c3(2)*Dp*exp(-lam2(4)*time);
p22 = p2_eq + c1(2)*App*exp(-lam2(2)*time) + c2(2)*Bpp*exp(-lam2(3)*time) + c3(2)*Dpp*exp(-lam2(4)*time);
p23 = p3_eq + c1(2)*Appp*exp(-lam2(2)*time) + c2(2)*Bppp*exp(-lam2(3)*time) + c3(2)*Dppp*exp(-lam2(4)*time);

p30 = p0_eq + c1(3)*A*exp(-lam2(2)*time) + c2(3)*B*exp(-lam2(3)*time) + c3(3)*D*exp(-lam2(4)*time);
p31 = p1_eq + c1(3)*Ap*exp(-lam2(2)*time) + c2(3)*Bp*exp(-lam2(3)*time) + c3(3)*Dp*exp(-lam2(4)*time);
p32 = p2_eq + c1(3)*App*exp(-lam2(2)*time) + c2(3)*Bpp*exp(-lam2(3)*time) + c3(3)*Dpp*exp(-lam2(4)*time);
p33 = p3_eq + c1(3)*Appp*exp(-lam2(2)*time) + c2(3)*Bppp*exp(-lam2(3)*time) + c3(3)*Dppp*exp(-lam2(4)*time);

tcf = A0*A0*p0_eq*p00 + A0*A1*p0_eq*p01 + A0*A2*p0_eq*p02 + A0*A3*p0_eq*p03;
tcf = tcf + A1*A0*p1_eq*p10 + A1*A1*p1_eq*p11 + A1*A2*p1_eq*p12 + A1*A3*p1_eq*p13;
tcf = tcf + A2*A0*p2_eq*p20 + A2*A1*p2_eq*p21 + A2*A2*p2_eq*p22 + A2*A3*p2_eq*p23;
tcf = tcf + A3*A0*p3_eq*p30 + A3*A1*p3_eq*p31 + A3*A2*p3_eq*p32 + A3*A3*p3_eq*p33;

%tcf = tcf/tcf(1);


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

end