function guess = knitro_gp32_3state_0p1uMgp32(x0)

obj = @(x) multigoaltcf_analytical(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));
% x(1), x(2), and x(3) are FRET values. x(4) through x(7) are fluctuation
% times.

% No contraints
nlcon = [];
AA = [];
bb = [];
Aeq = [];
beq = [];

% lower and upper bounds
lb = [.80999,.56,.55999,.0001,.001,.000001,.001,.001];
ub = [.81001,.8,.56001,.1,.1,.1,.1,.1];

%%%%%%%%%%%%%%%%%%% Parameters for computing rms %%%%%%%%%%%%%%%%%%%%%%%%%%
rmscalc = @multigoaltcf_analytical;

% Weights
w_2pt = 300000;
% w_4pt_t0 = 3000;
w_4pt_t1 = 5000;
% w_4pt_t3 = 2500;
w_4pt_t5 = 2500;
w_4pt_t20 = 2500;
w_4pt_t50 = 2500;
w_4pt_t100 = 2500;
w_frethist = 6000;%2000; % Weighting for FRET hist comparison % Original was 900. 9000 was crazy big.
w_eig1 = 3000;
w_eig2 = 3000;

tau1range = [1:250]';
tau3range = tau1range;

tcftime = [1:.1:500]';

% Weighting function
% w1func = ones(length(tcftime),1);
w1func = 1./sqrt(tcftime);
w2func = 1./sqrt(tau1range)*(1./sqrt(tau1range'));

% Open 0.1 uM gp32 experimental files
expt2pt = textread('TCFavgNorm_meansub_DoubleFitNew_0p1uMgp32.dat');
expt2pt = expt2pt(:,2);
expt2pt = expt2pt/expt2pt(1);
B = textread('FourCorrFinal_Fit_0p1uMgp32_tau1.dat');
B = B(tau1range,tau3range);
C = textread('FourCorrFinal_Fit_0p1uMgp32_tau5.dat');
C = C(tau1range,tau3range);
C1 = textread('FourCorrFinal_Fit_0p1uMgp32_tau20.dat');
C1 = C1(tau1range,tau3range);
C2 = textread('FourCorrFinal_Fit_0p1uMgp32_tau50.dat');
C2 = C2(tau1range,tau3range);
D = textread('FourCorrFinal_Fit_0p1uMgp32_tau100.dat');
D = D(tau1range,tau3range);

B = B/B(1,1);
C = C/C(1,1);
C1 = C1/C1(1,1);
C2 = C2/C2(1,1);
D = D/D(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End rms parameters %%%%%%%%%%%%%%%%%%%%

% solver call
[x,fit] = knitromatlab(obj, x0, AA, bb, Aeq, beq, lb, ub, nlcon);

guess = zeros(1,9);
guess(1:5) = x(4:8);
guess(6:8) = x(1:3);
guess(9) = fit;

function rms = multigoaltcf_analytical(val0,val1,val2,k01,k10,k20,k12,k21)
rms_array = zeros(15,1);

% Calculate 2-pt. TCF using analytical formula
tcf = TCF_cyclic3state_norm(tcftime,val0,val1,val2,k01,k10,k12,k21,k20);
tcf = tcf/tcf(1);

% Calculate rms deviation from 2-pt. TCF of guess
rms_array(1) = sum(((tcf - expt2pt).*w1func).^2)*w_2pt;

% Calculate 4-pt. TCF from analytical formula with various tau2 times.
[~,FourCorr1] = FourPtTCF_cyclic3state_norm(tau1range,1,val0,val1,val2,k01,k10,k12,k21,k20);
[~,FourCorr5] = FourPtTCF_cyclic3state_norm(tau1range,5,val0,val1,val2,k01,k10,k12,k21,k20);
[~,FourCorr20] = FourPtTCF_cyclic3state_norm(tau1range,20,val0,val1,val2,k01,k10,k12,k21,k20);
[~,FourCorr50] = FourPtTCF_cyclic3state_norm(tau1range,50,val0,val1,val2,k01,k10,k12,k21,k20);
[~,FourCorr100] = FourPtTCF_cyclic3state_norm(tau1range,100,val0,val1,val2,k01,k10,k12,k21,k20);

FourCorr1 = FourCorr1/FourCorr1(1,1);
FourCorr5 = FourCorr5/FourCorr5(1,1);
FourCorr20 = FourCorr20/FourCorr20(1,1);
FourCorr50 = FourCorr50/FourCorr50(1,1);
FourCorr100 = FourCorr100/FourCorr100(1,1);

% Calculate rms deviation from 4-pt. TCF of guess with t2 = 0
% rms_array(2) = sum(sum((A - FourCorr0).^2.*w2func))*w_4pt_t0;
rms_array(2) = sum(sum((B - FourCorr1).^2.*w2func))*w_4pt_t1;
% rms_array(3) = sum(sum((C - FourCorr3).^2.*w2func))*w_4pt_t3;
rms_array(4) = sum(sum((C - FourCorr5).^2.*w2func))*w_4pt_t5;
rms_array(5) = sum(sum((C1 - FourCorr20).^2.*w2func))*w_4pt_t20;
rms_array(6) = sum(sum((C2 - FourCorr50).^2.*w2func))*w_4pt_t50;
rms_array(7) = sum(sum((D - FourCorr100).^2.*w2func))*w_4pt_t100;
% rms_array(11) = sum(sum((C1 - FourCorr10).^2.*w2func))*w_4pt_t10;
% rms_array(12) = sum(sum((C2 - FourCorr15).^2.*w2func))*w_4pt_t15;
% rms_array(13) = sum(sum((D1 - FourCorr30).^2.*w2func))*w_4pt_t30;
% rms_array(14) = sum(sum((D2 - FourCorr40).^2.*w2func))*w_4pt_t40;
% rms_array(15) = sum(sum((E1 - FourCorr75).^2.*w2func))*w_4pt_t75;
rms_array(11) = 0;
rms_array(12) = 0;
rms_array(13) = 0;
rms_array(10) = 0;
rms_array(9) = 0;
rms_array(3) = 0;

% Comparing to 0.1 uM gp32 FRET histogram
bins = [0:.01:1]';
% fithist = 1336*exp(-((bins-.6528)/.08878).^2) + 1284*exp(-((bins-.663)/.1485).^2) + 348.4*exp(-((bins-.3181)/.09474).^2);
fithist = exp(-((bins-0.8045)/.09628).^2) + 753.8/2132*exp(-((bins-0.5682)/.1314).^2);
[p0_eq,p1_eq,p2_eq,eig1,eig2] = FourPtTCF_cyclic3state_xo(tau1range,0,val0,val1,val2,k01,k10,k12,k21,k20);
% val0scaled = 0.663;
% val3scaled = 0.318;
val0scaled = 0.8045;
val2scaled = 0.5682;
xscale = (val1 - val2)/(val0 - val2);
val1scaled = val2scaled + xscale*(val0scaled - val2scaled);
analhist = p0_eq*exp(-((bins-val0scaled)/.1).^2) + p1_eq*exp(-((bins-val1scaled)/.1).^2) + p2_eq*exp(-((bins-val2scaled)/.1).^2);
analhist = analhist/max(analhist);
rms_array(8) = sum(w_frethist*(fithist - analhist).^2);

% Compare to "known" eigenvalues
% Compare to "known" eigenvalues
if eig1 <= eig2
%     rms_array(14) = abs(eig2 - 0.1253)*w_eig1;
%     rms_array(15) = abs(eig1 - 0.01758)*w_eig2;
    rms_array(14) = abs(eig2 - 0.054347826086957)*w_eig1;
    rms_array(15) = abs(eig1 - 0.006337135614702)*w_eig2;
else
%     rms_array(14) = abs(eig1 - 0.1253)*w_eig1;
%     rms_array(15) = abs(eig2 - 0.01758)*w_eig2;
    rms_array(14) = abs(eig1 - 0.054347826086957)*w_eig1;
    rms_array(15) = abs(eig2 - 0.006337135614702)*w_eig2;
end

rms = sum(rms_array);

end

end