function [rms,rms_array] = Compare_4state_1p0uMgp32(guess)

% sv = 1;
sv = 0;
dir = 'C:\Users\Ultra-Carey\Desktop\Optimization Testing\BestGuess\FourState-1p0uM\JuneFitRoutine\Jun18th_knitro_Run1_histalt\';

a = 0.05848;
b = 0.8077;
c = 0.1157;
d = 0.5859;
e = 0.1895;
f = 368.1;

%val0 = guess(end,8);
%val1 = guess(end,9);
%val2 = guess(end,9);
%val3 = guess(end,10);
val0 = .81;
val3 = .56;
val1 = guess(end,9);
val2 = val1;
k01 = 1/guess(end,1);
k10 = 1/guess(end,2);
k20 = 1/guess(end,3);
k12 = 1/guess(end,4);
k21 = 1/guess(end,5);
k23 = 1/guess(end,6);
k32 = 1/guess(end,7);

k02 = k01*k12*k20/(k10*k21)

tau1range = [1:500]';
tau3range = [1:500]';

% Open 1.0 uM gp32 experimental files
% expt2pt = textread('TCFavg_1msRes_19Mol_meanSub.dat');
% expt2pt = expt2pt(:,2);
% expt2pt = expt2pt/(0.002701974600000);
expt2pt = textread('TCFavgNorm_meansub_DoubleFitNew_1p0uMgp32.dat');
expt2pt = expt2pt/expt2pt(1);
% expt2pt = expt2pt(1:4991);
% expt2ptfit = textread('TCFavgNorm_100usRes_1p0uM_DoubleFit.dat');
% expt2ptfit = expt2ptfit(:,2);
% A = textread('FourCorrFinal_1p0uMgp32_19mol_tau0.dat');
A = textread('FourCorrFinal_Fit_1p0uMgp32_tau0.dat');
A = A(tau1range,tau3range);
A0 = textread('FourCorrFinal_Fit_1p0uMgp32_tau1.dat');
A0 = A0(tau1range,tau3range);
B = textread('FourCorrFinal_Fit_1p0uMgp32_tau5.dat');
B = B(tau1range,tau3range);
B1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau10.dat');
B1 = B1(tau1range,tau3range);
B2 = textread('FourCorrFinal_Fit_1p0uMgp32_tau15.dat');
B2 = B2(tau1range,tau3range);
C = textread('FourCorrFinal_Fit_1p0uMgp32_tau20.dat');
C = C(tau1range,tau3range);
C1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau30.dat');
C1 = C1(tau1range,tau3range);
C2 = textread('FourCorrFinal_1p0uMgp32_19mol_tau40.dat');
C2 = C2(tau1range,tau3range);
D = textread('FourCorrFinal_1p0uMgp32_19mol_tau50.dat');
D = D(tau1range,tau3range);
D1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau75.dat');
D1 = D1(tau1range,tau3range);
E = textread('FourCorrFinal_1p0uMgp32_19mol_tau100.dat');
E = E(tau1range,tau3range);

% A = textread('FourCorrFinal_Fit_1p0uMgp32_tau0.dat');
% A = A(tau1range,tau3range);
% A0 = textread('FourCorrFinal_Fit_1p0uMgp32_tau1.dat');
% A0 = A0(tau1range,tau3range);
% B = textread('FourCorrFinal_Fit_1p0uMgp32_tau5.dat');
% B = B(tau1range,tau3range);
% B1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau10.dat');
% B1 = B1(tau1range,tau3range);
% B2 = textread('FourCorrFinal_Fit_1p0uMgp32_tau15.dat');
% B2 = B2(tau1range,tau3range);
% C = textread('FourCorrFinal_Fit_1p0uMgp32_tau20.dat');
% C = C(tau1range,tau3range);
% C1 = textread('FourCorrFinal_Fit_1p0uMgp32_tau30.dat');
% C1 = C1(tau1range,tau3range);
% C2 = textread('FourCorrFinal_1p0uMgp32_19mol_tau40.dat');
% C2 = C2(tau1range,tau3range);
% D = textread('FourCorrFinal_1p0uMgp32_19mol_tau50.dat');
% D = D(tau1range,tau3range);
% D1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau75.dat');
% D1 = D1(tau1range,tau3range);
% E = textread('FourCorrFinal_1p0uMgp32_19mol_tau100.dat');
% E = E(tau1range,tau3range);A = textread('FourCorrFinal_Fit_1p0uMgp32_tau0.dat');

A = A(tau1range,tau3range);
A0 = textread('FourCorrFinal_1p0uMgp32_19mol_tau1.dat');
A0 = A0(tau1range,tau3range);
B = textread('FourCorrFinal_1p0uMgp32_19mol_tau5.dat');
B = B(tau1range,tau3range);
B1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau10.dat');
B1 = B1(tau1range,tau3range);
B2 = textread('FourCorrFinal_1p0uMgp32_19mol_tau15.dat');
B2 = B2(tau1range,tau3range);
C = textread('FourCorrFinal_1p0uMgp32_19mol_tau20.dat');
C = C(tau1range,tau3range);
C1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau30.dat');
C1 = C1(tau1range,tau3range);
C2 = textread('FourCorrFinal_1p0uMgp32_19mol_tau40.dat');
C2 = C2(tau1range,tau3range);
D = textread('FourCorrFinal_1p0uMgp32_19mol_tau50.dat');
D = D(tau1range,tau3range);
D1 = textread('FourCorrFinal_1p0uMgp32_19mol_tau75.dat');
D1 = D1(tau1range,tau3range);
E = textread('FourCorrFinal_1p0uMgp32_19mol_tau100.dat');
E = E(tau1range,tau3range);

A = A/A(1,1);
A0 = A0/A0(1,1);
B = B/B(1,1);
B1 = B1/B1(1,1);
B2 = B2/B2(1,1);
C = C/C(1,1);
C1 = C1/C1(1,1);
C2 = C2/C2(1,1);
D = D/D(1,1);
D1 = D1/D1(1,1);
E = E/E(1,1);

% A = A(1:100,1:100);
% A0 = A0(1:100,1:100);
% B = B(1:100,1:100);
% B1 = B1(1:100,1:100);
% B2 = B2(1:100,1:100);
% C = C(1:100,1:100);
% C1 = C1(1:100,1:100);
% C2 = C2(1:100,1:100);
% D = D(1:100,1:100);
% D1 = D1(1:100,1:100);
% E = E(1:100,1:100);

% Calculate 2-pt. TCF from analytical formula.
tcftime = [1:.1:500]';
tcf = TCF_4state(tcftime,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
tcf = tcf/tcf(1);

% Calculate 4-pt. TCF from analytical formula with various tau2 times.
[~,FourCorr0] = FourPtTCF_FourState_v3(tau1range,0,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr1] = FourPtTCF_FourState_v3(tau1range,1,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr5] = FourPtTCF_FourState_v3(tau1range,5,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr10] = FourPtTCF_FourState_v3(tau1range,10,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr15] = FourPtTCF_FourState_v3(tau1range,15,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr20] = FourPtTCF_FourState_v3(tau1range,20,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr30] = FourPtTCF_FourState_v3(tau1range,30,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr40] = FourPtTCF_FourState_v3(tau1range,40,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr50] = FourPtTCF_FourState_v3(tau1range,50,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr75] = FourPtTCF_FourState_v3(tau1range,75,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
[~,FourCorr100] = FourPtTCF_FourState_v3(tau1range,100,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);

FourCorr0 = FourCorr0/FourCorr0(1,1);
FourCorr1 = FourCorr1/FourCorr1(1,1);
FourCorr5 = FourCorr5/FourCorr5(1,1);
FourCorr10 = FourCorr10/FourCorr10(1,1);
FourCorr15 = FourCorr15/FourCorr15(1,1);
FourCorr20 = FourCorr20/FourCorr20(1,1);
FourCorr30 = FourCorr30/FourCorr30(1,1);
FourCorr40 = FourCorr40/FourCorr40(1,1);
FourCorr50 = FourCorr50/FourCorr50(1,1);
FourCorr75 = FourCorr75/FourCorr75(1,1);
FourCorr100 = FourCorr100/FourCorr100(1,1);

figure
semilogx(tcftime,expt2pt,tcftime,tcf)
axis([0 300 -.1 1.1])
if sv == 1
    print('-dpng',[dir 'twopttcf'])
end

%%%%%%%%%%%%%%%%%%%%% New Plot Format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Parameters for computing rms %%%%%%%%%%%%%%%%%%%%%%%%%%
rmscalc = @multigoaltcf_analytical;

% Weights
w_2pt = 30000;
w_4pt_t0 = 0;
w_4pt_t1 = 1000;
w_4pt_t5 = 1000;
w_4pt_t10 = 1000;
w_4pt_t15 = 1000;
w_4pt_t20 = 1000;
w_4pt_t30 = 1000;
w_4pt_t40 = 0;%500;
w_4pt_t50 = 0;%500;
w_frethist = 10;

% Weighting functions
% 2-pt. TCF weighting function
w1func = 1./(tcftime);
% 4-pt. TCF weighting function
w2func = 1./tau1range*(1./(tau1range'));

[rms,rms_array] = rmscalc(val0,val1,val2,k01,k10,k20,k12,k21);
% rms = 0;
% rms_array = 0;

%%%%%%%%%%%%%%%%%%% End parameters for computing rms %%%%%%%%%%%%%%%%%%%%%%

%%%% 4-pt Image Plot %%%%%%%%%%%%
figure('Position', [100, 100, 1149, 295]);
subplot(2,9,1)
imagesc(FourCorr0)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,10)
imagesc(A)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,2)
imagesc(FourCorr1)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,11)
imagesc(A0)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,3)
imagesc(FourCorr5)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,12)
imagesc(B)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,4)
imagesc(FourCorr10)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,13)
imagesc(B1)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,5)
imagesc(FourCorr15)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,14)
imagesc(B2)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,6)
imagesc(FourCorr20)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,15)
imagesc(C)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,7)
imagesc(FourCorr30)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,16)
imagesc(C1)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,8)
imagesc(FourCorr40)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,17)
imagesc(C2)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,9)
imagesc(FourCorr50)
colormap jet
axis image
axis xy
axis([0 100 0 100])

subplot(2,9,18)
imagesc(D)
colormap jet
axis image
axis xy
axis([0 100 0 100])
if sv == 1
    print('-dpng',[dir 'FourCorrImages'])
end

%%%% End 4-Pt Image Plot %%%%%%%%%%%%

%%%% 4-pt Cross-section plots %%%%%%%%%%%%
figure
semilogx(tau1range,FourCorr0(:,1),'b',tau1range,FourCorr0(:,10),'r',tau1range,FourCorr0(:,50),'g',tau1range,FourCorr0(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 0 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,A(:,1),'b.',tau1range,A(:,10),'r.',tau1range,A(:,50),'g.',tau1range,A(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau0'])
end

figure
semilogx(tau1range,FourCorr1(:,1),'b',tau1range,FourCorr1(:,5),'m',tau1range,FourCorr1(:,10),'r',tau1range,FourCorr1(:,50),'g',tau1range,FourCorr1(:,100),'c',tau1range,FourCorr1(:,200),'b',tau1range,FourCorr1(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 1 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,A0(:,1),'b.',tau1range,A0(:,5),'m.',tau1range,A0(:,10),'r.',tau1range,A0(:,50),'g.',tau1range,A0(:,100),'c.',tau1range,A0(:,200),'b.',tau1range,A0(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau1'])
end

figure
semilogx(tau1range,FourCorr5(:,1),'b',tau1range,FourCorr5(:,5),'m',tau1range,FourCorr5(:,10),'r',tau1range,FourCorr5(:,50),'g',tau1range,FourCorr5(:,100),'c',tau1range,FourCorr5(:,200),'b',tau1range,FourCorr5(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 5 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,B(:,1),'b.',tau1range,B(:,5),'m.',tau1range,B(:,10),'r.',tau1range,B(:,50),'g.',tau1range,B(:,100),'c.',tau1range,B(:,200),'b.',tau1range,B(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau5'])
end

figure
semilogx(tau1range,FourCorr10(:,1),'b',tau1range,FourCorr10(:,5),'m',tau1range,FourCorr10(:,10),'r',tau1range,FourCorr10(:,50),'g',tau1range,FourCorr10(:,100),'c',tau1range,FourCorr10(:,200),'b',tau1range,FourCorr10(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 10 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,B1(:,1),'b.',tau1range,B1(:,5),'m.',tau1range,B1(:,10),'r.',tau1range,B1(:,50),'g.',tau1range,B1(:,100),'c.',tau1range,B1(:,200),'b.',tau1range,B1(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau10'])
end

figure
semilogx(tau1range,FourCorr15(:,1),'b',tau1range,FourCorr15(:,10),'r',tau1range,FourCorr15(:,50),'g',tau1range,FourCorr15(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 15 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,B2(:,1),'b.',tau1range,B2(:,10),'r.',tau1range,B2(:,50),'g.',tau1range,B2(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau15'])
end

figure
semilogx(tau1range,FourCorr20(:,1),'b',tau1range,FourCorr20(:,5),'m',tau1range,FourCorr20(:,10),'r',tau1range,FourCorr20(:,50),'g',tau1range,FourCorr20(:,100),'c',tau1range,FourCorr20(:,200),'b',tau1range,FourCorr20(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 20 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,C(:,1),'b.',tau1range,C(:,5),'m.',tau1range,C(:,10),'r.',tau1range,C(:,50),'g.',tau1range,C(:,100),'c.',tau1range,C(:,200),'b.',tau1range,C(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau20'])
end

figure
semilogx(tau1range,FourCorr30(:,1),'b',tau1range,FourCorr30(:,10),'r',tau1range,FourCorr30(:,50),'g',tau1range,FourCorr30(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 30 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,C1(:,1),'b.',tau1range,C1(:,10),'r.',tau1range,C1(:,50),'g.',tau1range,C1(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau30'])
end

figure
semilogx(tau1range,FourCorr40(:,1),'b',tau1range,FourCorr40(:,10),'r',tau1range,FourCorr40(:,50),'g',tau1range,FourCorr40(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 40 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,C2(:,1),'b.',tau1range,C2(:,10),'r.',tau1range,C2(:,50),'g.',tau1range,C2(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau40'])
end

figure
semilogx(tau1range,FourCorr50(:,1),'b',tau1range,FourCorr50(:,10),'r',tau1range,FourCorr50(:,50),'g',tau1range,FourCorr50(:,500),'k');
axis([0 500 -.05 1.05])
hold on
title('Compare fit to expt., tau2 = 50 ms')
xlabel('lag time (ms)')
ylabel('TCF')
semilogx(tau1range,D(:,1),'b.',tau1range,D(:,10),'r.',tau1range,D(:,50),'g.',tau1range,D(:,500),'k.')
if sv == 1
    print('-dpng',[dir 'FourCorr_tau50'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hst = textread('FREThist_1p0uMgp32.dat');
bins = [0:.01:1]';
[p0_eq,p1_eq,p2_eq,p3_eq] = FourPtTCF_FourState_xo(tau1range,1,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
val0scaled = 0.79;
val3scaled = 0.587;
xscale = (val1 - val3)/(val0 - val3);
val1scaled = val3scaled + xscale*(val0scaled - val3scaled);
val2scaled = val1scaled;
% analhist = p0_eq*exp(-((bins-val0scaled)/.1).^2) + p1_eq*exp(-((bins-val1scaled)/.1).^2) + p2_eq*exp(-((bins-val2scaled)/.1).^2) + p3_eq*exp(-((bins-val3scaled)/.1).^2);

%%% Best fit hist %%%
a = 0.05848;
b = 0.8077;
c = 0.1157;
d = 0.5859;
e = 0.1895;
f = 368.1;
analhist = f*(0.108476247855992/(a*sqrt(pi))*exp(-(bins-b).^2/(a^2)) + 0.543267825491906/(c*sqrt(pi))*exp(-(bins-d).^2/(c^2)) + 0.096062855687200/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)) + 0.252193070964902/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)));
%%%%%%%%%%%%%%%%%%%%%

% normsy = max(analhist);
% analhist = analhist/normsy;
figure
subplot(1,2,1)
bar(bins,analhist)
% axis([0 1 0 1.1]);
axis([0 1 0 1500]);
hold on
% plot(bins,p0_eq*exp(-((bins-val0scaled)/.1).^2)/normsy,bins,p1_eq*exp(-((bins-val1scaled)/.1).^2)/normsy,bins,p2_eq*exp(-((bins-val2scaled)/.1).^2)/normsy,bins,p3_eq*exp(-((bins-val3scaled)/.1).^2)/normsy,'LineWidth',3);
plot(bins,f*0.108476247855992/(a*sqrt(pi))*exp(-(bins-b).^2/(a^2)),bins,f*0.543267825491906/(c*sqrt(pi))*exp(-(bins-d).^2/(c^2)),bins,f*0.096062855687200/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)),bins,f*0.252193070964902/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)),'LineWidth',3);
subplot(1,2,2)
bar(bins,hst)
axis([0 1 0 1500]);
if sv == 1
    print('-dpng',[dir 'FREThist'])
end

% k02 = k01*k12*k20/(k10*k21);
% traj = FRETsim4state(30000,val0,val1,val2,val3,1/k01,1/k10,1/k02,1/k20,1/k12,1/k21,1/k23,1/k32);
% figure
% plot(traj)
% title('Simulated trajectory')
% xlabel('Time (ms)')
% ylabel('FRET')

function [rms,rms_array] = multigoaltcf_analytical(val0,val1,val2,k01,k10,k20,k12,k21)
rms_array = zeros(11,1);

% Calculate rms deviation from 2-pt. TCF of guess
rms_array(1) = sum(((tcf - expt2pt).*w1func).^2)*w_2pt;

% Calculate rms deviation from 4-pt. TCF of guess with t2 = 0
rms_array(2) = sum(sum((A - FourCorr0).^2.*w2func))*w_4pt_t0;
rms_array(3) = sum(sum((A0 - FourCorr1).^2.*w2func))*w_4pt_t1;
rms_array(4) = sum(sum((B - FourCorr5).^2.*w2func))*w_4pt_t5;
rms_array(5) = sum(sum((B1 - FourCorr10).^2.*w2func))*w_4pt_t10;
rms_array(6) = sum(sum((B2 - FourCorr15).^2.*w2func))*w_4pt_t15;
rms_array(7) = sum(sum((C - FourCorr20).^2.*w2func))*w_4pt_t20;
rms_array(8) = sum(sum((C1 - FourCorr30).^2.*w2func))*w_4pt_t30;
rms_array(9) = sum(sum((C2 - FourCorr40).^2.*w2func))*w_4pt_t40;
rms_array(10) = sum(sum((C2 - FourCorr50).^2.*w2func))*w_4pt_t50;

% Comparing to 1.0 uM FRET histogram
bins = [0:.01:1]';
fithist = 686.9/1271*exp(-((bins-0.8045)/.09628).^2) + exp(-((bins-0.5682)/.1314).^2);
[p0_eq,p1_eq,p2_eq,p3_eq,~,~,~] = FourPtTCF_FourState_xo(tau1range,1,val0,val1,val2,val3,k01,k10,k20,k12,k21,k23,k32);
val0scaled = 0.8045;
val2scaled = 0.5682;
xscale = (val1 - val2)/(val0 - val2);
val1scaled = val2scaled + xscale*(val0scaled - val2scaled);
% analhist = f*(0.108476247855992/(a*sqrt(pi))*exp(-(bins-b).^2/(a^2)) + 0.543267825491906/(c*sqrt(pi))*exp(-(bins-d).^2/(c^2)) + 0.096062855687200/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)) + 0.252193070964902/(e*sqrt(pi))*exp(-(bins-d-.472*(b-d)).^2/(e^2)));
analhist = p0_eq*exp(-((bins-val0scaled)/.1).^2) + p1_eq*exp(-((bins-val1scaled)/.1).^2) + p2_eq*exp(-((bins-val1scaled)/.1).^2) + p3_eq*exp(-((bins-val2scaled)/.1).^2);
analhist = analhist/max(analhist);
rms_array(11) = sum(w_frethist*(fithist - analhist).^2);

rms = sum(rms_array);

end

end