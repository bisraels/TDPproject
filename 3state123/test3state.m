%--------------------------------------------------------------------------
% User Options
%--------------------------------------------------------------------------

clockMode = 1;
plotMode = 1;
verboseMode = 0;

plotCondProbMode = 1;
plot_histMode = 1;
plot_C2Mode = 1;
plot_C4Mode = 1;
%--------------------------------------------------------------------------
% Build the symbolic conditional probabilities
%--------------------------------------------------------------------------
%Number of states
N = 3;
if verboseMode == 1
    disp(['N = ' num2str(N)]);
end

% Define symbolic matrix of c's for our constants
c = sym('c', [N,N]);
%c =
% [ c1_1, c1_2, c1_3]
% [ c2_1, c2_2, c2_3]
% [ c3_1, c3_2, c3_3]
if verboseMode == 1
    disp('c = ');
    disp(c);
end

% Define symbolic matrix of the eigenvector components.
% Transpose the matrix for proPCr matrix mulitplication
v = sym('v',[N,N]).';
if verboseMode == 1
    disp('v = ');
    disp(v);
end
% v =
% [ v1_1, v2_1, v3_1]
% [ v1_2, v2_2, v3_2]
% [ v1_3, v2_3, v3_3]

% Make a symbolic vector of the eigenvalues and time
lam = exp(sym('eval',[N,1])*'time');

if verboseMode == 1
    disp('lam = ');
    disp(lam);
end
% lam_sym =
%  exp(lam1*time)
%  exp(lam2*time)
%  exp(lam3*time)

%make a symbolic figure
p = v*(c.*lam);
if verboseMode == 1
    pline = p(:);
    disp('pline = ')
    disp(pline);
end


% pline =
%  c1_1*v1_1*exp(lam1*time) + c2_1*v2_1*exp(lam2*time) + c3_1*v3_1*exp(lam3*time)
%  c1_1*v1_2*exp(lam1*time) + c2_1*v2_2*exp(lam2*time) + c3_1*v3_2*exp(lam3*time)
%  c1_1*v1_3*exp(lam1*time) + c2_1*v2_3*exp(lam2*time) + c3_1*v3_3*exp(lam3*time)
%  c1_2*v1_1*exp(lam1*time) + c2_2*v2_1*exp(lam2*time) + c3_2*v3_1*exp(lam3*time)
%  c1_2*v1_2*exp(lam1*time) + c2_2*v2_2*exp(lam2*time) + c3_2*v3_2*exp(lam3*time)
%  c1_2*v1_3*exp(lam1*time) + c2_2*v2_3*exp(lam2*time) + c3_2*v3_3*exp(lam3*time)
%  c1_3*v1_1*exp(lam1*time) + c2_3*v2_1*exp(lam2*time) + c3_3*v3_1*exp(lam3*time)
%  c1_3*v1_2*exp(lam1*time) + c2_3*v2_2*exp(lam2*time) + c3_3*v3_2*exp(lam3*time)
%  c1_3*v1_3*exp(lam1*time) + c2_3*v2_3*exp(lam2*time) + c3_3*v3_3*exp(lam3*time)

%--------------------------------------------------------------------------
% Pick random values for the rate constants and FRET states
%--------------------------------------------------------------------------
[k12,k13,k21,k23,k31,A1,A2,A3,k32,~] = paramSim_3state123_cyclical();

% k12 = 12; k13 = 13; k21 = 21; k31 = 31; k23 = 23;
% k32 = k12*k23*k31/(k13*k21);
% A1 = .1; A2 = .2; A3 = .3;

% load('BestFitResults_fmincon.mat','A1','A2','A3','k12','k13','k21','k23','k31','k32')
%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors of K matrix
%--------------------------------------------------------------------------
K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];

[Vec, lambda] = eig(K);

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

% Assign numerical values to the N eigenvalue symbols
eval1 = double(lambdaSorted(1,1));
eval2 = double(lambdaSorted(2,2));
eval3 = double(lambdaSorted(3,3));
evals = sum(lambdaSorted);
%lam = sum(lambdaSorted);

evec1 = VecSorted(:,1);
evec2 = VecSorted(:,2);
evec3 = VecSorted(:,3);


v1_1 = evec1(1);
v1_2 = evec1(2);
v1_3 = evec1(3);

v2_1 = evec2(1);
v2_2 = evec2(2);
v2_3 = evec2(3);

v3_1 = evec3(1);
v3_2 = evec3(2);
v3_3 = evec3(3);

% Npts = 150;
% time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;


%Load the algebraic solution for the expansion coefficients
load('wCoef_3state_condensed_sol.mat',...
    'c1_1', 'c1_2', 'c1_3',...
    'c2_1', 'c2_2', 'c2_3',...
    'c3_1', 'c3_2', 'c3_3');

%Substitute in the eigenvector coefficents into the expresssions for the
%expansion coefficients
c1_1 = double(subs(c1_1));
c2_1 = double(subs(c2_1));
c3_1 = double(subs(c3_1));

c1_2 = double(subs(c1_2));
c2_2 = double(subs(c2_2));
c3_2 = double(subs(c3_2));

c1_3 = double(subs(c1_3));
c2_3 = double(subs(c2_3));
c3_3 = double(subs(c3_3));



%At time t = infinity
P1_eq = c1_1*v1_1;
P2_eq = c1_2*v1_2;
P3_eq = c1_3*v1_3;


PCq = [P1_eq, P2_eq, P3_eq];
eqSum = sum(PCq);
disp(['The sum of equilibrium probabilities is ' num2str(eqSum)]);

A = [A1,A2,A3];



time = time_sim();
%Set the conditional probabilities using the full expressions
% % Pj_i = c1_i*v1_j*exp(lam1*itme) + c2_i*v2_j*exp(lam2*itme) + c3_i*v3_j*exp(lam3*itme);
p1_1 = c1_1*v1_1*exp(eval1*time) + c2_1*v2_1*exp(eval2*time) + c3_1*v3_1*exp(eval3*time);
p2_1 = c1_1*v1_2*exp(eval1*time) + c2_1*v2_2*exp(eval2*time) + c3_1*v3_2*exp(eval3*time);
p3_1 = c1_1*v1_3*exp(eval1*time) + c2_1*v2_3*exp(eval2*time) + c3_1*v3_3*exp(eval3*time);
p1_2 = c1_2*v1_1*exp(eval1*time) + c2_2*v2_1*exp(eval2*time) + c3_2*v3_1*exp(eval3*time);
p2_2 = c1_2*v1_2*exp(eval1*time) + c2_2*v2_2*exp(eval2*time) + c3_2*v3_2*exp(eval3*time);
p3_2 = c1_2*v1_3*exp(eval1*time) + c2_2*v2_3*exp(eval2*time) + c3_2*v3_3*exp(eval3*time);
p1_3 = c1_3*v1_1*exp(eval1*time) + c2_3*v2_1*exp(eval2*time) + c3_3*v3_1*exp(eval3*time);
p2_3 = c1_3*v1_2*exp(eval1*time) + c2_3*v2_2*exp(eval2*time) + c3_3*v3_2*exp(eval3*time);
p3_3 = c1_3*v1_3*exp(eval1*time) + c2_3*v2_3*exp(eval2*time) + c3_3*v3_3*exp(eval3*time);
% % note: pj_eq = c1_i*v1_j (i.e. p2_eq = c1_1*v1_2 = c1_2*v1_2 = c1_3*v1_2)


%-------------------------------------------------------------------------
% Calculate the conditional probabilities using matrices
%-------------------------------------------------------------------------

clear('time');

%Simulate a sym array
lam = exp(sym('eval',[N,1])*'time');

%Fill in lam with eigenvalues variables (NOT TIME)
Lam1 = subs(lam(1));
Lam2 = subs(lam(2));
Lam3 = subs(lam(3));
LamVec = [Lam1;Lam2;Lam3];

C = [...
    c1_1, c1_2, c1_3;...
    c2_1, c2_2, c2_3;...
    c3_1, c3_2, c3_3;...
    ];
if verboseMode == 1
    disp('C = ');
    disp(C);
end

V = VecSorted;
% V = [...
%     v1_1, v2_1, v3_1;...
%     v1_2, v2_2, v3_2;...
%     v1_3, v2_3, v3_3;
%     ];
if verboseMode == 1
    disp('V = ');
    disp(V);
end

%Make a matrix of conditional probabilites as a funtion of time
p_t_sym = V*(C.*LamVec);

%transform the NxN matrix into a (N*N)x1 matrix
p_t_sym_line = p_t_sym(:);

%Define a time
time = time_sim();

%Replace pline_sym with the time
Pline = double(subs(p_t_sym_line));

%Recast 9xnumel(time) matrix back into original size (NxNxnumel(t))
P = reshape(Pline, [N N numel(time)]);


if plot_histMode == 1
    figure(5);
    clf
    
    FRET_bins = [0:.01:1];
    sigma_A1 = 0.15;
    sigma_A2 = 0.1;
    sigma_A3 = 0.1;
    
    
    Peq = [P(1,1,end), P(2,2,end), P(3,3,end)];
    p1_eq = Peq(1);
    p2_eq = Peq(2);
    p3_eq = Peq(3);
    
    hist_sim = p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2) + p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2) + p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2);
    denom_hist_sim = sum(hist_sim);
    hist_sim = hist_sim./sum(hist_sim);
    hold on;
    
    fitHistPlot = plot(FRET_bins,hist_sim,'DisplayName','Cumulative Fit');
    fitHistPlot.LineStyle = '-';
    %         fitHistPlot.Color = 'red';
    fitHistPlot.LineWidth = 2;
    
    state1_HistPlot = plot(FRET_bins,p1_eq*exp(-((FRET_bins-A1)/sigma_A1).^2)./denom_hist_sim,...
        'k--','LineWidth',1,'DisplayName',['p1_{eq} =' num2str(p1_eq)]);
    state2_HistPlot = plot(FRET_bins,p2_eq*exp(-((FRET_bins-A2)/sigma_A2).^2)./denom_hist_sim,...
        'b--','LineWidth',1,'DisplayName',['p2_{eq} =' num2str(p2_eq)]);
    state3_HistPlot = plot(FRET_bins,p3_eq*exp(-((FRET_bins-A3)/sigma_A3).^2)./denom_hist_sim,...
        'r--','LineWidth',1,'DisplayName',['p3_{eq} =' num2str(p3_eq)]);
   
    
    
        xlabel('FRET Efficiency','FontSize',14);
        ylabel('Frequency','FontSize',14);
        
        [sample_description, ~] = sample_descriptionGetter();
        title_str = ['Experimental vs Simulated Histograms' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        
        lgd = legend('show');
        lgd.Location = 'northwest';
        lgd.FontSize = 14;
        drawnow();
        hold on;
        whitebg;
        set(gca,'FontSize',14);
end

%%
%-------------------------------------------------------------------------
% Plot the Conditional Probabilities
%-------------------------------------------------------------------------
if plotCondProbMode == 1
    if plotMode == 1
        
        figure(3);
        clf;
        set(gcf,'Color','w');
        hold on;
        logx;
        set(gca,'xscale','log');
        
        [sample_description, ~] = sample_descriptionGetter();
        title_str = ['Conditional Proababilities: P_{i\rightarrowj}(t)' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        %         title('Conditional Proababilities','FontSize',18)
        xlabel('Time (sec)','FontSize',14);
        ylabel('Probability','FontSize',14);
        
        % 1 := C    Coiled
        % 2 := PC   Partially Coiled
        % 3 := FE   Fully Extended
        
        state_1_color = 'k';%[0, 0.4470, 0.7410];        % "C"
        state_2_color = 'b';%[0.8500, 0.3250, 0.0980];   % "I"
        state_3_color = 'r';%[0.9290, 0.6940, 0.1250];   % "FE"
        
        timeUB = 1e-2;
        %Begin in state 1 -block
        plot(time,p1_1,'color',state_1_color,'LineWidth',3,'LineStyle','-','DisplayName','P_{{\color{black}C}}\rightarrow_{C}');
        text(timeUB*.8,p1_1(end)+0.05,'P_C^{Eq}','FontSize',18);
        plot(time,p2_1,'color',state_1_color,'LineWidth',3,'LineStyle',':','DisplayName','P_{C}\rightarrow_{I}');
        plot(time,p3_1,'color',state_1_color,'LineWidth',3,'LineStyle','--','DisplayName','P_{C}\rightarrow_{FE}')
        
          lgd = legend;
        lgd.FontSize = 14;
        lgd.Location = 'bestoutside';
        xlim([-inf,timeUB]);
        set(gca,'FontSize',14);
        
        
        %Begin in state-2
        plot(time,p1_2,'color',state_2_color,'LineWidth',3,'LineStyle','-','DisplayName','P_{I}\rightarrow_{C}');
        plot(time,p2_2,'color',state_2_color,'LineWidth',3,'LineStyle',':','DisplayName','P_{I}\rightarrow_{I}');
        text(timeUB*.8,p2_2(end)+0.05,'P_I^{Eq}','FontSize',18);
        plot(time,p3_2,'color',state_2_color,'LineWidth',3,'LineStyle','--','DisplayName','P_{I}\rightarrow_{FE}');
        
        
        %Begin in state-3
        plot(time,p1_3,'color',state_3_color,'LineWidth',3,'LineStyle','-','DisplayName','P_{FE}\rightarrow_{C}');
        plot(time,p2_3,'color',state_3_color,'LineWidth',3,'LineStyle',':','DisplayName','P_{FE}\rightarrow_{I}');
        plot(time,p3_3,'color',state_3_color,'LineWidth',3,'LineStyle','--','DisplayName','P_{FE}\rightarrow_{FE}');
        text(timeUB*.8,p3_3(end)+0.05,'P_{FE}^{Eq}','FontSize',18);
        
        
      
    end
end
%%
if plot_C2Mode == 1
%--------------------------------------------------------------------------
% Calculate the Two-point TCF
%--------------------------------------------------------------------------
% Calculate C2_sim (METHOD 1)
% tic
% C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
%     P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
%     P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
% disp('     Time to calculate C2 using one line...');
% toc

% Calculate C2_sim (METHOD 2)
%mean has already been subtracted
A = [A1,A2,A3];
Amean = sum(A.*PCq);
% Subtract mean values
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

Peq = [P1_eq,P2_eq,P3_eq];
%
% cP = zeros(numel(A),numel(A),length(time));
% cP(1,1,:) = p1_1;
% cP(2,1,:) = p2_1;
% cP(3,1,:) = p3_1;
% cP(1,3,:) = p1_2;
% cP(2,3,:) = p2_2;
% cP(3,3,:) = p3_2;
% cP(1,3,:) = p1_3;
% cP(2,3,:) = p2_3;
% cP(3,3,:) = p3_3;

% tic
% C2_sim2 = zeros(size(time));
% for i = 1:numel(A)
%     for j = 1:numel(A)
%         C2_sim2_temp = A(j) * squeeze(cP(j,i,:)) * A(i) * PCq(i);
%
%         %         C2_sim1_temp = A(j) * reshape(cP(j,i,:),numel(cP)) * A(i) * PCq(i);
%         C2_sim2 = C2_sim2 + C2_sim2_temp;
%     end
% end
% disp('     Time to calculate C2 using loops and squeeze...');
% toc

% % Calculate C2_sim (METHOD 3)
% tic
% C2_sim3 = zeros(size(time));
% for i = 1:numel(A)
%     for j = 1:numel(A)
%         C2_sim3_temp = A(j) * squeeze(P(j,i,:)) * A(i) * PCq(i);
%         %         C2_sim1_temp = A(j) * reshape(cP(j,i,:),numel(cP)) * A(i) * PCq(i);
%         C2_sim3 = C2_sim3 + C2_sim3_temp;
%     end
% end
% disp('     Time to calculate C2 using loops and squeeze (method 3)...');
% toc

% Calculate C2_sim (METHOD 4)
tic
Npts = numel(time);
C2_sim3 = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        C2_sim3_temp = A(j) * reshape(P(j,i,:),[1,Npts]) * A(i) * PCq(i);
        %         C2_sim1_temp = A(j) * reshape(cP(j,i,:),numel(cP)) * A(i) * PCq(i);
        C2_sim3 = C2_sim3 + C2_sim3_temp;
    end
end
disp('     Time to calculate C2 using loops and reshape (method 4)...');
toc
%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------
if plotMode == 1
    figure(2)
    clf;
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,3,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for sPCcfic time range
    %     TCF2pt = plot(time,C2_sim,'LineWidth',3,'DisplayName','C2 one line');
    hold on;
    %     TCF2pt2 = plot(time,C2_sim2,'r--','LineWidth',3,'DisplayName','C2 loops');
    TCF2pt3 = plot(time,C2_sim3,'b-','LineWidth',3,'DisplayName','C2 loops matrix');
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    
    ax = gca;
    ax.XScale = 'log';
    
end
end

%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------
if plot_C4Mode == 1
% Calculate the NxN conditional Probabilities for tau2
tau2 = 0;
p1_1_t2 = P1_eq + c2_1*v2_1*exp(eval2*tau2) + c3_1*v3_1*exp(eval3*tau2);
p1_2_t2 = P1_eq + c2_2*v2_1*exp(eval2*tau2) + c3_2*v3_1*exp(eval3*tau2);
p1_3_t2 = P1_eq + c2_3*v2_1*exp(eval2*tau2) + c3_3*v3_1*exp(eval3*tau2);
p2_1_t2 = P2_eq + c2_1*v2_2*exp(eval2*tau2) + c3_1*v3_2*exp(eval3*tau2);
p2_2_t2 = P2_eq + c2_2*v2_2*exp(eval2*tau2) + c3_2*v3_2*exp(eval3*tau2);
p2_3_t2 = P2_eq + c2_3*v2_2*exp(eval2*tau2) + c3_3*v3_2*exp(eval3*tau2);
p3_1_t2 = P3_eq + c2_1*v2_3*exp(eval2*tau2) + c3_1*v3_3*exp(eval3*tau2);
p3_2_t2 = P3_eq + c2_2*v2_3*exp(eval2*tau2) + c3_2*v3_3*exp(eval3*tau2);
p3_3_t2 = P3_eq + c2_3*v2_3*exp(eval2*tau2) + c3_3*v3_3*exp(eval3*tau2);

cP_t2 = zeros(numel(A),numel(A),length(time));
cP_t2(1,1,:) = p1_1_t2 ;
cP_t2(2,1,:) = p2_1_t2 ;
cP_t2(3,1,:) = p3_1_t2 ;
cP_t2(1,3,:) = p1_2_t2 ;
cP_t2(2,3,:) = p2_2_t2 ;
cP_t2(3,3,:) = p3_2_t2 ;
cP_t2(1,3,:) = p1_3_t2 ;
cP_t2(2,3,:) = p2_3_t2 ;
cP_t2(3,3,:) = p3_3_t2 ;

%-------------------------------------------------------------------------
% Iterate over all the Permutations of FRET States
%-------------------------------------------------------------------------
% tic
% C4 = zeros(numel(time),numel(time));
% for i = 1:numel(A)
%     for j = 1:numel(A)
%         for k = 1:numel(A)
%             for l = 1:numel(A)
%                 C4term_val =  A(l) *squeeze(cP(l,k,:)) * A(k) * cP_t2(k,j) * A(j) * squeeze(cP(j,i,:))'* A(i) * PCq(i);
%                 C4 = C4 + C4term_val;
%             end
%         end
%     end
% end
% disp('    time to calculate C4 using loops and squeeze');
% toc

%-------------------------------------------------------------------------
% Iterate over all the PCrmutations of FRET States (MAtrix Methods)
%-------------------------------------------------------------------------
%Define a time
time_holder = time;
time = 0;

%Replace pline_sym with the time
Pline_tau2 = double(subs(p_t_sym_line));
%Recast back into original size
P_tau2 = reshape(Pline_tau2, [N N numel(time)]);
time = time_holder;
% 
% tic
% C4_sim2 = zeros(numel(time),numel(time));
% for i = 1:numel(A)
%     for j = 1:numel(A)
%         for k = 1:numel(A)
%             for l = 1:numel(A)
%                 C4term_val =  A(l) *squeeze(P(l,k,:)) * A(k) * P_tau2(k,j) * A(j) * squeeze(P(j,i,:))'* A(i) * PCq(i);
%                 C4_sim2 = C4_sim2 + C4term_val;
%             end
%         end
%     end
% end
% disp('    time to calculate C4 using loops and squeeze (Matrix Method)');
% toc
% 

tic
Npts = numel(time);
C4_sim3 = zeros(numel(time),numel(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                C4term_val =  A(l) *reshape(P(l,k,:),[1,Npts])' * A(k) *  P_tau2(k,j) * A(j) * reshape(P(j,i,:),[1,numel(time)])* A(i) * PCq(i);
                C4_sim3 = C4_sim3 + C4term_val;
            end
        end
    end
end
disp('    time to calculate C4 using loops and reshape (Matrix Method)');
toc

%-------------------------------------------------------------------------
% Plot the 4 point TCF
%-------------------------------------------------------------------------
if plotMode == 1
    
    figure(4);
    clf;
%     surf(time, time, C4,'DisplayName','Conditional Prob');
    hold on
    
%     mesh(time,time,C4_sim2,'DisplayName','loops/squeeze');
    mesh(time,time,C4_sim3,'DisplayName','loops reshape');
    
    legend();
    title('Four-point TCF: C^{(4)}','FontSize',18)
    xlabel('Time (\tau_1)','FontSize',14);
    ylabel('Time (\tau_3)','FontSize',14);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
    
    set(gcf,'color','w');
    
    view(28,36);
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
    drawnow();
    hold on;
end
end
