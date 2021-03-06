%exp: Can we use N^2 eqns? Yes but its slower. 
%Can we eliminate all but one of the expansion coefficients? No.
% Can we solve 3 equations at once? All Pji for a certain j


clockMode = 1;
plotMode = 1;
verboseMode = 1;
% Define matrix of c's for our constants
c = sym('c', [3,3]);

% Define matrix of the eigenvector components
V = sym('v',[3,3]).' ;

% Use identity matrix to define the possible conditions
cond = eye(3);

% Use these to construct our equations
eq1 = cond(:,1) == V * c(:,1) ;
eq2 = cond(:,2) == V * c(:,2) ;
eq3 = cond(:,3) == V * c(:,3) ;

% Define the eqilibrium probabilities
P1_eq = (c(1,1) * V(1,1)); % from P11
P2_eq = (c(1,2) * V(2,1)); % from P22
P3_eq = (c(1,3) * V(3,1)); % from P33


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
% % % %    tic
% % % %    
% % % % c1_1 = solsc1.c1_1;
% % % % c2_1 = solsc1.c2_1;
% % % % c3_1 = solsc1.c3_1;
% % % % 
% % % % toc
% % % % 
% % % % tic
% % % % c1_2 = solsc2.c1_2;
% % % % c2_2 = solsc2.c2_2;
% % % % c3_2 = solsc2.c3_2;
% % % % 
% % % % toc
% % % % tic
% % % % c1_3 = solsc3.c1_3;
% % % % c2_3 = solsc3.c2_3;
% % % % c3_3 = solsc3.c3_3;
% % % %    toc
% % % % 


%%
% Solve eqns for vars
tic
sols = vpasolve(eqs, vars);
disp(['     Time to solve ' num2str(numel(eqs)) ' coupled equations ...']);
toc
%Takes about 0.3 seconds


% Plug in the solutions to solve for the eqilibrium probabilities
P1_eq = (sols.c1_1 * V(1,1)); 
P2_eq = (sols.c1_2 * V(2,1)); 
P3_eq = (sols.c1_3 * V(3,1)); 

c1_1 = sols.c1_1;
c2_1 = sols.c2_1;
c3_1 = sols.c3_1;

c1_2 = sols.c1_2;
c2_2 = sols.c2_2;
c3_2 = sols.c3_2;

c1_3 = sols.c1_3;
c2_3 = sols.c2_3;
c3_3 = sols.c3_3;

% save('wCoef_3state_condensed.mat','sols',...
%     'c1_1', 'c1_2', 'c1_3',...
%     'c2_1', 'c2_2', 'c2_3',...
%     'c3_1', 'c3_2', 'c3_3',...
%     'P1_eq', 'P2_eq','P3_eq');


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


% if clockMode == 1
%     tic
% end
[Vec, lambda] = eig(K);   
% %Display the amount of time a process took. Begins at the last tic.
% if clockMode == 1
%     elapsedTime = toc;
%     task_str = 'Calculate the eigenvalues';
%     disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
% end

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
% % % 
% % % P1_eq = double(subs(P1_eq));
% % % P2_eq = double(subs(P2_eq));
% % % P3_eq = double(subs(P3_eq));
% % % 
% % % % Calculate the NxN conditional Probabilities
% % % p1_1 = P1_eq + c2_1*v2_1*exp(lam2*time) + c3_1*v3_1*exp(lam3*time);
% % % p1_2 = P1_eq + c2_2*v2_1*exp(lam2*time) + c3_2*v3_1*exp(lam3*time);
% % % p1_3 = P1_eq + c2_3*v2_1*exp(lam2*time) + c3_3*v3_1*exp(lam3*time);
% % % p2_1 = P2_eq + c2_1*v2_2*exp(lam2*time) + c3_1*v3_2*exp(lam3*time);
% % % p2_2 = P2_eq + c2_2*v2_2*exp(lam2*time) + c3_2*v3_2*exp(lam3*time);
% % % p2_3 = P2_eq + c2_3*v2_2*exp(lam2*time) + c3_3*v3_2*exp(lam3*time);
% % % p3_1 = P3_eq + c2_1*v2_3*exp(lam2*time) + c3_1*v3_3*exp(lam3*time);
% % % p3_2 = P3_eq + c2_2*v2_3*exp(lam2*time) + c3_2*v3_3*exp(lam3*time);
% % % p3_3 = P3_eq + c2_3*v2_3*exp(lam2*time) + c3_3*v3_3*exp(lam3*time);
% % % 
% % % % Compute the 2 point TCF (Using conditinal Probabilities)
% % % % Subtract mean values
% % % Amean = P1_eq*A1 + P2_eq*A2 + P3_eq*A3;
% % % A1 = A1 - Amean;
% % % A2 = A2 - Amean;
% % % A3 = A3 - Amean;
% % % 
% % % 
% % % 
% % % % Calculate C2_sim
% % % tic
% % % C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
% % %     P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
% % %     P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
% % % disp('     Time to calculate C2 using one line...');
% % % toc
% % % % Matrix of conditional probabilities Pji (i-->j) with i is initial condition
% % % % cP = [ p1_1, p1_2, p1_3;
% % % %     p2_1, p2_2, p2_3;
% % % %     p3_1, p3_2, p3_3];
% % % 
% % % %mean has already been subtracted
% % % A = [A1,A2,A3];
% % % Peq = [P1_eq,P2_eq,P3_eq];
% % % 
% % % cP = zeros(numel(A),numel(A),length(time));
% % % cP(1,1,:) = p1_1;
% % % cP(2,1,:) = p2_1;
% % % cP(3,1,:) = p3_1;
% % % cP(1,2,:) = p1_2;
% % % cP(2,2,:) = p2_2;
% % % cP(3,2,:) = p3_2;
% % % cP(1,3,:) = p1_3;
% % % cP(2,3,:) = p2_3;
% % % cP(3,3,:) = p3_3;
% % % 
% % % %--------------------------------------------------------------------------
% % % tic
% % % C2_sim2 = zeros(size(time));
% % % for i = 1:numel(A)
% % %     for j = 1:numel(A)
% % %                 C2_sim1_temp = A(j) * squeeze(cP(j,i,:)) * A(i) * Peq(i);
% % % 
% % % %         C2_sim1_temp = A(j) * reshape(cP(j,i,:),numel(cP)) * A(i) * Peq(i);
% % %         C2_sim2 = C2_sim2 + C2_sim1_temp;
% % %     end
% % % end
% % % disp('     Time to calculate C2 using loops and squeeze...');
% % % toc
% % % 
% % % %--------------------------------------------------------------------------
% % % %  Plot two point TCF
% % % %--------------------------------------------------------------------------
% % % if plotMode == 1
% % %     figure(2)
% % %     clf;
% % %     set(gcf,'Color','w');
% % %     set(gcf,'Name','C2');
% % %     %subplot(1,2,1)
% % %     %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
% % %     TCF2pt = plot(time,C2_sim,'LineWidth',2);
% % %     hold on;
% % %     TCF2pt2 = plot(time,C2_sim2,'r--','LineWidth',2);
% % %     
% % %     title('Two point TCF','FontSize',18)
% % %     xlabel('Time','FontSize',14);
% % %     ylabel('C^{(2)}(\tau)','FontSize',14);
% % % 
% % %     ax = gca;
% % %     ax.XScale = 'log';
% % %     
% % % end
% % % 
% % % % Calculate the NxN conditional Probabilities
% % % tau2 = 0;
% % % p1_1_t2 = P1_eq + c2_1*v2_1*exp(lam2*tau2) + c3_1*v3_1*exp(lam3*tau2);
% % % p1_2_t2 = P1_eq + c2_2*v2_1*exp(lam2*tau2) + c3_2*v3_1*exp(lam3*tau2);
% % % p1_3_t2 = P1_eq + c2_3*v2_1*exp(lam2*tau2) + c3_3*v3_1*exp(lam3*tau2);
% % % p2_1_t2 = P2_eq + c2_1*v2_2*exp(lam2*tau2) + c3_1*v3_2*exp(lam3*tau2);
% % % p2_2_t2 = P2_eq + c2_2*v2_2*exp(lam2*tau2) + c3_2*v3_2*exp(lam3*tau2);
% % % p2_3_t2 = P2_eq + c2_3*v2_2*exp(lam2*tau2) + c3_3*v3_2*exp(lam3*tau2);
% % % p3_1_t2 = P3_eq + c2_1*v2_3*exp(lam2*tau2) + c3_1*v3_3*exp(lam3*tau2);
% % % p3_2_t2 = P3_eq + c2_2*v2_3*exp(lam2*tau2) + c3_2*v3_3*exp(lam3*tau2);
% % % p3_3_t2 = P3_eq + c2_3*v2_3*exp(lam2*tau2) + c3_3*v3_3*exp(lam3*tau2);
% % % 
% % % cP = zeros(numel(A),numel(A),length(time));
% % % cP_t2(1,1,:) = p1_1_t2 ;
% % % cP_t2(2,1,:) = p2_1_t2 ;
% % % cP_t2(3,1,:) = p3_1_t2 ;
% % % cP_t2(1,2,:) = p1_2_t2 ;
% % % cP_t2(2,2,:) = p2_2_t2 ;
% % % cP_t2(3,2,:) = p3_2_t2 ;
% % % cP_t2(1,3,:) = p1_3_t2 ;
% % % cP_t2(2,3,:) = p2_3_t2 ;
% % % cP_t2(3,3,:) = p3_3_t2 ;
% % % 
% % % %-------------------------------------------------------------------------
% % % % Iterate over all the Permutations of FRET States
% % % %-------------------------------------------------------------------------
% % % tic
% % % C4 = zeros(numel(time),numel(time));
% % % for i = 1:numel(A)
% % %     for j = 1:numel(A)
% % %         for k = 1:numel(A)
% % %             for l = 1:numel(A)
% % %                 C4term_val =  A(l) *squeeze(cP(l,k,:)) * A(k) * cP_t2(k,j) * A(j) * squeeze(cP(j,i,:))'* A(i) * Peq(i);
% % %                 C4 = C4 + C4term_val;
% % %             end
% % %         end
% % %     end
% % % end
% % % disp('    time to calculate C4');
% % % toc

%-------------------------------------------------------------------------
% Calculate the conditional probabilities using matrices
%-------------------------------------------------------------------------

%Make a matrix C of expansion coefficient components
C =[...
    c1_1, c1_2, c1_3;...
   c2_1, c2_2, c2_3;...
   c3_1, c3_2, c3_3;...
   ];
if verboseMode == 1
disp('C = ');
    disp(C)
end
V = VecSorted;

N = 6;
clear('time');

lam = exp(sym('eval',[1,N])*'time');
Lam1 =  subs(lam(1));
Lam2 =  subs(lam(2));
Lam3 =  subs(lam(3));
Lam4 =  subs(lam(4));
Lam5 =  subs(lam(5));
Lam6 =  subs(lam(6));
LamVec = [Lam1;Lam2;Lam3;Lam4;Lam5;Lam6];

%Make a matrix of conditional probabilites as a funtion of time
p_notime_sym = V*(C.*LamVec);

%transform the NxN matrix into a (N*N)x1 matrix
p_notime_sym_line = p_notime_sym(:);

%Define a time
time = time_sim();

%Replace pline_sym with the time
Pline = double(subs(p_notime_sym_line));

%Recast 9xnumel(time) matrix back into original size (NxNxnumel(t))
P = reshape(Pline, [N N numel(time)]);

%--------------------------------------------------------------------------
% Calculate C2 using the 3x3xnumel(time) sized conditional probability
% Matrix P where P(j,i,:) is the prob of going from i->j in time t
%--------------------------------------------------------------------------
tic
Npts = numel(time);
C2_sim = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        C2_sim_temp = A(j) * reshape(P(j,i,:),[1,Npts]) * A(i) * Peq(i);
        C2_sim = C2_sim + C2_sim_temp;
    end
end
elapsedTime = toc;
disp(['Took ' num2str(elapsedTime) ' seconds to calculate C2 for N = ' num2str(N)]);
%Takes about 2 msec

%--------------------------------------------------------------------------
%  Plot two point TCF
%--------------------------------------------------------------------------
if plotMode == 1
    figure(2)
    clf;
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    TCF2pt = plot(time,C2_sim,'LineWidth',2,'DisplayName','C2');
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);
    
    ax = gca;
    ax.XScale = 'log';
    
end


%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------

tic
Npts = numel(time);
C4_sim = zeros(numel(time),numel(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                C4term_val =  A(l) *reshape(P(l,k,:),[1,Npts])' * A(k) *  P_tau2(k,j) * A(j) * reshape(P(j,i,:),[1,numel(time)])* A(i) * Peq(i);
                C4_sim = C4_sim + C4term_val;
            end
        end
    end
end
elapsedTime = toc;
disp(['Took ' num2str(elapsedTime) ' seconds to calculate C4 for N = ' num2str(N)]);
%Takes about 7 msec

%-------------------------------------------------------------------------
% Plot the 4 point TCF
%-------------------------------------------------------------------------
if plotMode == 1
    figure(4);
    clf;
    surf(time, time, C4_sim,'DisplayName','Conditional Prob Matrix');
    hold on
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