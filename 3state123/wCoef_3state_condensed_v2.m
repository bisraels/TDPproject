
%--------------------------------------------------------------------------
% User Options
%--------------------------------------------------------------------------


clockMode = 1;
plotMode = 1;
saveMode = 0;

% Define matrix of c's for our constants
c = sym('c', [3,3]);
%c = 
% [ c1_1, c1_2, c1_3]
% [ c2_1, c2_2, c2_3]
% [ c3_1, c3_2, c3_3]

% Define matrix of the eigenvector components
v = sym('v',[3,3]).' ;
% v = 
% [ v1_1, v2_1, v3_1]
% [ v1_2, v2_2, v3_2]
% [ v1_3, v2_3, v3_3]

% Use identity matrix to define the possible conditions
cond = eye(3);

% Use these to construct our equations
eq1 = cond(:,1) == v * c(:,1) ;
eq2 = cond(:,2) == v * c(:,2) ;
eq3 = cond(:,3) == v * c(:,3) ;

% Define the eqilibrium probabilities
P1_eq = (c(1,1) * v(1,1)); % from P11
P2_eq = (c(1,2) * v(2,1)); % from P22
P3_eq = (c(1,3) * v(3,1)); % from P33


% Additional constraint such that Peq sums to 1.
eq_sum = P1_eq + P2_eq + P3_eq == 1;

% eqs = [eq1(1),eq1(2),eq1(3),...
%     eq2(1),eq2(2),eq2(3),...
%     eq3(1),eq3(2),eq3(3),...
%     eq_sum];
eqs = [eq1, eq2, eq3, eq_sum];

% Build list of vars - writing out by hand is prone to  mistakes
vars = c(:);

% if length(vars) == (length(c)*length(c))
%     disp('We have all the variables!')
% else
%     disp('We have missing/extra variables.')
% end
   

% Solve eqns for vars
tic
sols = vpasolve(eqs, vars);
disp(['     Time to solve ' num2str(numel(eqs)) ' coupled equations ...']);
toc
%Takes about 0.3 seconds


% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * v(1,1)); 
P2_eq = (sols.c1_2 * v(2,1)); 
P3_eq = (sols.c1_3 * v(3,1)); 

c1_1 = sols.c1_1;
c2_1 = sols.c2_1;
c3_1 = sols.c3_1;

c1_2 = sols.c1_2;
c2_2 = sols.c2_2;
c3_2 = sols.c3_2;

c1_3 = sols.c1_3;
c2_3 = sols.c2_3;
c3_3 = sols.c3_3;

if saveMode == 1
    fnameOut = 'wCoef_3state_condensed_sol.mat';
    save(fnameOut,'sols',...
        'c1_1', 'c1_2', 'c1_3',...
        'c2_1', 'c2_2', 'c2_3',...
        'c3_1', 'c3_2', 'c3_3',...
        'P1_eq', 'P2_eq','P3_eq');
    disp(['Saved the solution as ' fnameOut]);
end

load('wCoef_3state_condensed_sol.mat',...
    'c1_1', 'c1_2', 'c1_3',...
    'c2_1', 'c2_2', 'c2_3',...
    'c3_1', 'c3_2', 'c3_3',...
    'P1_eq', 'P2_eq','P3_eq');
% Test the solutions with calcuation of 2 point TCF
%--------------------------------------------------------------------------
% Pick random values for the rate constants and FRET states
%--------------------------------------------------------------------------
[k12,k13,k21,k23,k31,A1,A2,A3,k32,time] = paramSim_3state123_cyclical();

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors
%--------------------------------------------------------------------------
K = [(-k12 - k13), k21, k31;...
    k12, (-k21 - k23 ), k32;...
    k13, k23, (-k31-k32);];


if clockMode == 1
    tic
end
[Vec, lambda] = eig(K);   
%Display the amount of time a process took. Begins at the last tic.
if clockMode == 1
    elapsedTime = toc;
    task_str = 'Calculate the eigenvalues';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

lam1 = double(lambdaSorted(1,1));
lam2 = double(lambdaSorted(2,2));
lam3 = double(lambdaSorted(3,3));

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

% Calculate the expansion coefficients in terms of the eigevcetor
% components

c1_1 = double(subs(c1_1));
c2_1 = double(subs(c2_1));
c3_1 = double(subs(c3_1));

c1_2 = double(subs(c1_2));
c2_2 = double(subs(c2_2));
c3_2 = double(subs(c3_2));

c1_3 = double(subs(c1_3));
c2_3 = double(subs(c2_3));
c3_3 = double(subs(c3_3));

P1_eq = double(subs(P1_eq));
P2_eq = double(subs(P2_eq));
P3_eq = double(subs(P3_eq));

% Calculate the NxN conditional Probabilities
p1_1 = P1_eq + c2_1*v2_1*exp(lam2*time) + c3_1*v3_1*exp(lam3*time);
p1_2 = P1_eq + c2_2*v2_1*exp(lam2*time) + c3_2*v3_1*exp(lam3*time);
p1_3 = P1_eq + c2_3*v2_1*exp(lam2*time) + c3_3*v3_1*exp(lam3*time);
p2_1 = P2_eq + c2_1*v2_2*exp(lam2*time) + c3_1*v3_2*exp(lam3*time);
p2_2 = P2_eq + c2_2*v2_2*exp(lam2*time) + c3_2*v3_2*exp(lam3*time);
p2_3 = P2_eq + c2_3*v2_2*exp(lam2*time) + c3_3*v3_2*exp(lam3*time);
p3_1 = P3_eq + c2_1*v2_3*exp(lam2*time) + c3_1*v3_3*exp(lam3*time);
p3_2 = P3_eq + c2_2*v2_3*exp(lam2*time) + c3_2*v3_3*exp(lam3*time);
p3_3 = P3_eq + c2_3*v2_3*exp(lam2*time) + c3_3*v3_3*exp(lam3*time);

% Compute the 2 point TCF (Using conditinal Probabilities)
% Subtract mean values
Amean = P1_eq*A1 + P2_eq*A2 + P3_eq*A3;
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;



% % Calculate C2_sim
% tic
% C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
%     P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
%     P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
% disp('     Time to calculate C2 using one line...');
% toc

%mean has already been subtracted
A = [A1,A2,A3];
Peq = [P1_eq,P2_eq,P3_eq];

% Matrix of conditional probabilities Pji (i-->j) with i is initial condition
% cP = [ p1_1, p1_2, p1_3;
%     p2_1, p2_2, p2_3;
%     p3_1, p3_2, p3_3];
cP_time = zeros(numel(A),numel(A),length(time));
cP_time(1,1,:) = p1_1;
cP_time(2,1,:) = p2_1;
cP_time(3,1,:) = p3_1;
cP_time(1,2,:) = p1_2;
cP_time(2,2,:) = p2_2;
cP_time(3,2,:) = p3_2;
cP_time(1,3,:) = p1_3;
cP_time(2,3,:) = p2_3;
cP_time(3,3,:) = p3_3;

%--------------------------------------------------------------------------
tic
C2_sim2 = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
                C2_sim1_temp = A(j) * squeeze(cP_time(j,i,:)) * A(i) * Peq(i);

%         C2_sim1_temp = A(j) * reshape(cP(j,i,:),numel(cP)) * A(i) * Peq(i);
        C2_sim2 = C2_sim2 + C2_sim1_temp;
    end
end
disp('     Time to calculate C2 using loops and squeeze...');
toc

%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------
if plotMode == 1
    figure(2)
    clf;
    set(gcf,'position',[170 501 560 420]);
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,2,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
%     TCF2pt = plot(time,C2_sim,'LineWidth',2);
    hold on;
    TCF2pt2 = plot(time,C2_sim2,'b','LineWidth',2);
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);

    ax = gca;
    ax.XScale = 'log';
    
end

%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


% Calculate the NxN conditional Probabilities for tau2
tau2 = 0;
p1_1_t2 = P1_eq + c2_1*v2_1*exp(lam2*tau2) + c3_1*v3_1*exp(lam3*tau2);
p1_2_t2 = P1_eq + c2_2*v2_1*exp(lam2*tau2) + c3_2*v3_1*exp(lam3*tau2);
p1_3_t2 = P1_eq + c2_3*v2_1*exp(lam2*tau2) + c3_3*v3_1*exp(lam3*tau2);
p2_1_t2 = P2_eq + c2_1*v2_2*exp(lam2*tau2) + c3_1*v3_2*exp(lam3*tau2);
p2_2_t2 = P2_eq + c2_2*v2_2*exp(lam2*tau2) + c3_2*v3_2*exp(lam3*tau2);
p2_3_t2 = P2_eq + c2_3*v2_2*exp(lam2*tau2) + c3_3*v3_2*exp(lam3*tau2);
p3_1_t2 = P3_eq + c2_1*v2_3*exp(lam2*tau2) + c3_1*v3_3*exp(lam3*tau2);
p3_2_t2 = P3_eq + c2_2*v2_3*exp(lam2*tau2) + c3_2*v3_3*exp(lam3*tau2);
p3_3_t2 = P3_eq + c2_3*v2_3*exp(lam2*tau2) + c3_3*v3_3*exp(lam3*tau2);

cP_t2_time = zeros(numel(A),numel(A),length(time));
cP_t2_time(1,1,:) = p1_1_t2 ;
cP_t2_time(2,1,:) = p2_1_t2 ;
cP_t2_time(3,1,:) = p3_1_t2 ;
cP_t2_time(1,2,:) = p1_2_t2 ;
cP_t2_time(2,2,:) = p2_2_t2 ;
cP_t2_time(3,2,:) = p3_2_t2 ;
cP_t2_time(1,3,:) = p1_3_t2 ;
cP_t2_time(2,3,:) = p2_3_t2 ;
cP_t2_time(3,3,:) = p3_3_t2 ;

%-------------------------------------------------------------------------
% Iterate over all the Permutations of FRET States
%-------------------------------------------------------------------------
tic
C4 = zeros(numel(time),numel(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                C4term_val =  A(l) *squeeze(cP_time(l,k,:)) * A(k) * cP_t2_time(k,j) * A(j) * squeeze(cP_time(j,i,:))'* A(i) * Peq(i);
                C4 = C4 + C4term_val;
            end
        end
    end
end
disp('    time to calculate C4 using loops and squeeze');
toc

%-------------------------------------------------------------------------
 % Plot the 4 point TCF
 %-------------------------------------------------------------------------
 if plotMode == 1
        figure(4);
        clf
         set(gcf,'position',[845 452 512 467]);
        surf(time, time, C4);
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        
        drawnow();
        hold on;
 end
    