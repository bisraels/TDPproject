clockMode = 1;
plotMode = 1;


%Number of states
N = 3;

% Define symbolic matrix of c's for our constants
c = sym('c', [N,N])
%c = 
% [ c1_1, c1_2, c1_3]
% [ c2_1, c2_2, c2_3]
% [ c3_1, c3_2, c3_3]

% Define symbolic matrix of the eigenvector components.
% Transpose the matrix for proper matrix mulitplication
v = sym('v',[N,N]).'
% v = 
% [ v1_1, v2_1, v3_1]
% [ v1_2, v2_2, v3_2]
% [ v1_3, v2_3, v3_3]

% Make a symbolic vector of the eigenvalues and time
lam_sym = exp(sym('lam',[N,1])*'time')
% lam_sym =
%  exp(lam1*time)
%  exp(lam2*time)
%  exp(lam3*time)

p_sym = v*(c.*lam_sym);
pline_sym = p_sym(:)

% pline_sym = 
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
[k12,k13,k21,k23,k31,A1,A2,A3,k32,time] = paramSim_3state123_cyclical();

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
lam1 = double(lambdaSorted(1,1));
lam2 = double(lambdaSorted(2,2));
lam3 = double(lambdaSorted(3,3));
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

%Fill in lam with eigenvalues and time variables
Lam1 = double(subs(lam_sym(1)));
Lam2 = double(subs(lam_sym(2)));
Lam3 = double(subs(lam_sym(3)));
LamVec = [Lam1;Lam2;Lam3]


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


 C =[c1_1, c1_2, c1_3;...
   c2_1, c2_2, c2_3, ;...
   c3_1, c3_2, c3_3, ;]
disp(C)

V = VecSorted.';
% time = 1:10;
% P = V*(C.*LamVec)
% lam_sym
P_notime_sym = V*(C.*lam_sym)


% P = P(:)

% substitute in numerical values for expansion coeff. and eigenvectors and
% time
% time = 0:11;
% pline = p(:)

%Replace pline_sym with the time
Pline = double(subs(P_notime_sym));
%Recast back into original size
P = reshape(Pline, [N N numel(time)]);

%At time t = infinity
P1_eq = c1_1*v1_1;
P2_eq = c1_2*v1_2;
P3_eq = c1_3*v1_3;


Peq = [P1_eq, P2_eq, P3_eq];
eqSum = sum(Peq);
disp(['The sum of equilibrium probabilities is ' num2str(eqSum)]);

A = [A1,A2,A3];
Amean = sum(A.*Peq);

% Subtract mean values
A1 = A1 - Amean;
A2 = A2 - Amean;
A3 = A3 - Amean;

%Set the conditional probabilities using the full expressions
% % Pj_i = c1_i*v1_j*exp(lam1*itme) + c2_i*v2_j*exp(lam2*itme) + c3_i*v3_j*exp(lam3*itme);
  p1_1 = c1_1*v1_1*exp(lam1*time) + c2_1*v2_1*exp(lam2*time) + c3_1*v3_1*exp(lam3*time);
  p2_1 = c1_1*v1_2*exp(lam1*time) + c2_1*v2_2*exp(lam2*time) + c3_1*v3_2*exp(lam3*time);
  p3_1 = c1_1*v1_3*exp(lam1*time) + c2_1*v2_3*exp(lam2*time) + c3_1*v3_3*exp(lam3*time);
  p1_2 = c1_2*v1_1*exp(lam1*time) + c2_2*v2_1*exp(lam2*time) + c3_2*v3_1*exp(lam3*time);
  p2_2 = c1_2*v1_2*exp(lam1*time) + c2_2*v2_2*exp(lam2*time) + c3_2*v3_2*exp(lam3*time);
  p3_2 = c1_2*v1_3*exp(lam1*time) + c2_2*v2_3*exp(lam2*time) + c3_2*v3_3*exp(lam3*time);
  p1_3 = c1_3*v1_1*exp(lam1*time) + c2_3*v2_1*exp(lam2*time) + c3_3*v3_1*exp(lam3*time);
  p2_3 = c1_3*v1_2*exp(lam1*time) + c2_3*v2_2*exp(lam2*time) + c3_3*v3_2*exp(lam3*time);
  p3_3 = c1_3*v1_3*exp(lam1*time) + c2_3*v2_3*exp(lam2*time) + c3_3*v3_3*exp(lam3*time);
% % note: pj_eq = c1_i*v1_j (i.e. p2_eq = c1_1*v1_2 = c1_2*v1_2 = c1_3*v1_2)


% Calculate C2_sim
tic
C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
    P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
    P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);
disp('     Time to calculate C2 using one line...');
toc


%mean has already been subtracted
A = [A1,A2,A3];
Peq = [P1_eq,P2_eq,P3_eq];

cP = zeros(numel(A),numel(A),length(time));
cP(1,1,:) = p1_1;
cP(2,1,:) = p2_1;
cP(3,1,:) = p3_1;
cP(1,2,:) = p1_2;
cP(2,2,:) = p2_2;
cP(3,2,:) = p3_2;
cP(1,3,:) = p1_3;
cP(2,3,:) = p2_3;
cP(3,3,:) = p3_3;

%--------------------------------------------------------------------------
tic
C2_sim2 = zeros(size(time));
for i = 1:numel(A)
    for j = 1:numel(A)
                C2_sim1_temp = A(j) * squeeze(cP(j,i,:)) * A(i) * Peq(i);

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
    set(gcf,'Color','w');
    set(gcf,'Name','C2');
    %subplot(1,2,1)
    %TCF2pt = fplot(C2(t),[1e-3,1],'LineWidth',2);      % fplot() was making it hard to plot on loglog scale, so calculate for specfic time range
    TCF2pt = plot(time,C2_sim,'LineWidth',2);
    hold on;
    TCF2pt2 = plot(time,C2_sim2,'r--','LineWidth',2);
    
    title('Two point TCF','FontSize',18)
    xlabel('Time','FontSize',14);
    ylabel('C^{(2)}(\tau)','FontSize',14);

    ax = gca;
    ax.XScale = 'log';
    
end


%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------

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

cP_t2 = zeros(numel(A),numel(A),length(time));
cP_t2(1,1,:) = p1_1_t2 ;
cP_t2(2,1,:) = p2_1_t2 ;
cP_t2(3,1,:) = p3_1_t2 ;
cP_t2(1,2,:) = p1_2_t2 ;
cP_t2(2,2,:) = p2_2_t2 ;
cP_t2(3,2,:) = p3_2_t2 ;
cP_t2(1,3,:) = p1_3_t2 ;
cP_t2(2,3,:) = p2_3_t2 ;
cP_t2(3,3,:) = p3_3_t2 ;

%-------------------------------------------------------------------------
% Iterate over all the Permutations of FRET States
%-------------------------------------------------------------------------
tic
C4 = zeros(numel(time),numel(time));
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                C4term_val =  A(l) *squeeze(cP(l,k,:)) * A(k) * cP_t2(k,j) * A(j) * squeeze(cP(j,i,:))'* A(i) * Peq(i);
                C4 = C4 + C4term_val;
            end
        end
    end
end
disp('    time to calculate C4');
toc

%-------------------------------------------------------------------------
 % Plot the 4 point TCF
 %-------------------------------------------------------------------------
 if plotMode == 1
        figure(4);
        surf(time, time, C4);
        title('Four-point TCF: C^{(4)}','FontSize',18)
        xlabel('Time (\tau_1)','FontSize',14);
        ylabel('Time (\tau_3)','FontSize',14);
        zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','FontSize',14);
        
        whitebg;
        
        view(28,36);
        ax = gca;
        ax.XScale = 'log';
        ax.YScale = 'log';
        
        drawnow();
        hold on;
 end