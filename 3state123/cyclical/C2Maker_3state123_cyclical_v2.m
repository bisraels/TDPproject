%--------------------------------------------------------------------------
% AUTHOR: Claire Albrecht & Brett Israels
%
% CREATED: September 2019 (C2Maker_3state123_cyclical.m)
%
% PURPOSE:  Evaluate the cyclical 3-State (123) conditional probabilties with a set of
%           rates then the  2 point TCF and 4 point TCF for that system
%
% INPUT: (1) conditional probabilities from ODE solver: symCondProb_3state123_cyclical.mat
%
% OUTPUT: (1) FRET Histogram
%         (2) Two-point TCF C2
%         (3) Four-point TCF C4
%
% MODIFICATIONS:
%   (1) Adapted from function_3state123_cyclical.m
%
%--------------------------------------------------------------------------

% function C2_sim = C2Maker_3state123_cyclical(t12,t13,t21,t23,t31,A1,A2,A3,time)
%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
verboseMode = 0; %Set to 1 to see alot of progress updates and print off.
clockMode = 1;
saveMode = 0;
plotMode = 0;

programName = 'C2Maker_3state123_cyclical';
%--------------------------------------------------------------------------
% SET PARAMATERS
%--------------------------------------------------------------------------
% switch nargin
%     case 0
%--------------------------------------------------------------------------
% Pick random values for the rate constants and FRET states
%--------------------------------------------------------------------------
[k12,k13,k21,k23,k31,A1,A2,A3,k32,time] = paramSim_3state123_cyclical();

        
%     case 8
%         
%         Npts = 150;
%         time = [0:9,logspace(1,6.4771212,Npts)]/1e6;
%         
%         
% end
%--------------------------------------------------------------------------
% Set the rates(5 rates)
%--------------------------------------------------------------------------
% k12 = 1/t12;
% k13 = 1/t13;
% k21 = 1/t21;
% k23 = 1/t23;
% k31 = 1/t31;
% 
% % % Detailed balance condition: %k31 will be the rate fixed by the others
% k32 = k12*k23*k31/(k13*k21);
% t32 = 1/k32;

%--------------------------------------------------------------------------
% Define the FRET Array
%--------------------------------------------------------------------------
A = [ A1; A2; A3];

%--------------------------------------------------------------------------
% User Prefrences
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
%Output created by ODE solver
if verboseMode == 1
    disp('Loading the conditional Probabilities as a function of rates');
end
load('symCondProb_3state123_cyclical.mat','P11','P12','P13','P21','P22','P23','P31','P32','P33','eval1','eval2','eval3')

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'load the conditional probabilities Pij';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

if clockMode == 1
    tic
end
% disp('Calculating conditional probabilities using the rates defined')
t = sym('t');

%--------------------------------------------------------------------------
% CALCULATE THE EIGENVALUES
%--------------------------------------------------------------------------

if clockMode == 1
    tic
end
% Evaluate the eigenvalues in terms of the rates defined above - produce as doubles
%subs(s) returns a copy of s, replacing symbolic variables in s, with their
%values obtained from the calling function and the MATLAB® Workspace,
% and then evaluates s. Variables with no assigned values remain as variables.
eval1 = double(vpa(subs(eval1)));
eval2 = double(vpa(subs(eval2)));
eval3 = double(vpa(subs(eval3)));

if verboseMode == 1
    disp(['The final eigentimescales are: tau1 = 1/lam1 = ' num2str(1e6*1/eval1) ' microseconds '...
        ' tau2 = 1/eval2 = ' num2str(1e6*1/eval2) ' microseconds ' ...
        'and tau3 = 1/eval3 = ' num2str(1e6*1/eval3) ' microseconds .']);
end

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'evaluate the eigenvalues using subs';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end


%--------------------------------------------------------------------------
% Evaluate the conditional probabilities
%--------------------------------------------------------------------------

if clockMode == 1
    tic
end
% Evaluate conditional probabilties by substituting in values from above
% and using vpa() to force the simplest form of the output.

%Pij(t) is prob from i--> j, assuming you start in state i: Pi(t=0)=100%=1

% P11(t) = vpa(subs(P11));
% P12(t) = vpa(subs(P12));
% P13(t) = vpa(subs(P13));
%
% P21(t) = vpa(subs(P21));
% P22(t) = vpa(subs(P22));
% P23(t) = vpa(subs(P23));
%
% P31(t) = vpa(subs(P31));
% P32(t) = vpa(subs(P32));
% P33(t) = vpa(subs(P33));

P11(t) = subs(P11);
P12(t) = subs(P12);
P13(t) = subs(P13);

P21(t) = subs(P21);
P22(t) = subs(P22);
P23(t) = subs(P23);

P31(t) = subs(P31);
P32(t) = subs(P32);
P33(t) = subs(P33);

%Slowed it down signifigantly
% P11 = subs(P11);
% P12 = subs(P12);
% P13 = subs(P13);
% 
% P21 = subs(P21);
% P22 = subs(P22);
% P23 = subs(P23);
% 
% P31 = subs(P31);
% P32 = subs(P32);
% P33 = subs(P33);


%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'evaluate the conditional probabilities as a function of rates {kij} using subs';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end
%
%--------------------------------------------------------------------------
% CALCULATE Equilibrium populations
%--------------------------------------------------------------------------

% t = sym('t');   % This allows for a dynamic workspace

t = time;
% Matrix of conditional probabilities Pi-->j with i is initial condition
cP = [ P11(t), P21(t), P31(t);
    P12(t), P22(t), P32(t);
    P13(t), P23(t), P33(t)];
% Row = final state & Column = initial condition
% Define equilibrium populations from the conditional probabiltiies at
% infinite time.
P1EQ = P11(inf);
P2EQ = P22(inf);
P3EQ = P33(inf);

Peq = [P1EQ; P2EQ; P3EQ];

%--------------------------------------------------------------------------
% (2) Calculate Two point TCF:
%-------------------------------------------------------------------------
Amean = sum(A.*Peq);
A = A - Amean;

if clockMode == 1
    tic
end

% C2sym(t) = 0*t;
C2sim = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        %
%         C2temp(t) = A(j) * cP(j,i) * A(i) * Peq(i);
%         C2sym(t) = C2sym(t) + C2temp(t);
C2temp  = A(j) * cP(j,i) * A(i) * Peq(i);
C2sim = C2sim + C2temp;
    end
end

%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the 2-point TCF as a function of rates {kij}';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

msq = sum((A.^2).*Peq);     % square of mean <A^2>.
sqm = (sum(A.*Peq))^2;      % mean square value <A>^2


if verboseMode == 1
    if double(C2sym(0)) == double(msq)
        disp('Mean of the square <A^2> matches C2(t=0)!')
    else
        fprintf('Problem: mean of the square <A^2> (%f) DOES NOT match C2(t=0) (%f)',double(msq),double(C2sym(0)))
    end
    
    if double(C2sym(10^20)) == double(sqm)
        disp('Square of the mean <A>^2 matches C2(t=inf)!')
    else
        
        fprintf('Problem: square of the mean <A>^2 (%f) DOES NOT match C2(t=inf) (%f)',double(sqm),double(C2sym(10^20)))
    end
end

%--------------------------------------------------------------------------
% Evaluate C2 over a range of t's
%--------------------------------------------------------------------------
% C2_sim = C2sym(time);

%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------

%close all
if plotMode == 1
    figure(2)
    
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,2,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
    TCF2pt = plot(time,C2_sim,'LineWidth',2);
    
    title('Analytical Two point TCF','FontSize',18)
    xlabel('Time (\tau_1)','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    %xlim([10^-5 10^0])
    % ylim([C2(500) C2(0)])
    
    ax = gca;
    ax.XScale = 'log';
    set(gca,'yscale','log')
    
    if saveMode == 1
        saveName = ['C2_','example'];
        saveas(TCF2pt,saveName, 'png')
    end
end

