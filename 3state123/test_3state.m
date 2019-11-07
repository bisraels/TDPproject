clockMode = 1;
plotMode = 1;


%Number of states
N = 3

% Define symbolic matrix of c's for our constants
c = sym('c', [N,N])

% Define symbolic matrix of the eigenvector components
v = sym('v',[N,N]).' 

% Make a symbolic vector of the eigenvalues and time
lam = exp(sym('lam',[1,N])*'time')

% time = 1:10;
p = v*(c.*lam.')


%Recast back into original size
% P2 = reshape(Pline, Psize)


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


%Load the algebraic solution for the expansion coefficients
load('wCoef_3state_condensed.mat',...
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




% P = P(:)

% substitute in numerical values for expansion coeff. and eigenvectors and
% time
time = 0:11;
pline = p(:)
Pline_numbers = double(subs(pline));
%Recast back into original size
CP = reshape(Pline_numbers, [N N numel(time)]);

%At time t = infinity
P1_eq = c1_1*v1_1
P2_eq = c1_2*v1_2
P3_eq = c1_3*v1_3


P = [P1_eq, P2_eq, P3_eq];
eqSum = sum(P);
disp(['The sum of equilibrium probabilities is ' num2str(eqSum)]);

A = [A1,A2,A3];
Amean = sum(A.*P);

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