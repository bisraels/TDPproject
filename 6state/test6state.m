
%--------------------------------------------------------------------------
% User Options
%--------------------------------------------------------------------------
clockMode = 1;
plotMode = 1;
verboseMode = 0;

plotCondProbMode = 0;


%--------------------------------------------------------------------------
% Build the symbolic conditional probabilities
%--------------------------------------------------------------------------
%Number of states
N = 6;
if verboseMode == 1
    disp(['N = ' num2str(N)]);
end

% Define symbolic matrix of c's for our constants
c = sym('c', [N,N]);
if verboseMode == 1
    disp('c = ');
    disp(c);
end

% Define symbolic matrix of the eigenvector components.
% Transpose the matrix for proper matrix mulitplication
v = sym('v',[N,N]).';
if verboseMode == 1
    disp('v = ');
    disp(v);
end
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


%--------------------------------------------------------------------------
% Simulate a K matrix for a 6-state system
%--------------------------------------------------------------------------
K = [-25,21,31,0,0,0;...
    12,-93,2852/91,42,52,0;...
    13,23,-(31+2852/91),0,0,0;...
    0,24,0,-42,0,0;...
    0,25,0,0,-108,65;...
    0,0,0,0,56,-65];

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors of K matrix
%--------------------------------------------------------------------------
tic

[Vec, lambda] = eig(K);

elapsedTime = toc;
if clockMode == 1
    task_str = ['solve the ' num2str(N) 'N-state eigensystem.'];
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);

eval1 = double(lambdaSorted(1,1));
eval2 = double(lambdaSorted(2,2));
eval3 = double(lambdaSorted(3,3));
eval4 = double(lambdaSorted(4,4));
eval5 = double(lambdaSorted(5,5));
eval6 = double(lambdaSorted(6,6));

evec1 = VecSorted(:,1);
evec2 = VecSorted(:,2);
evec3 = VecSorted(:,3);
evec4 = VecSorted(:,4);
evec5 = VecSorted(:,5);
evec6 = VecSorted(:,6);

%No getting around assigning variables one number at a time
v1_1 = evec1(1);
v1_2 = evec1(2);
v1_3 = evec1(3);
v1_4 = evec1(4);
v1_5 = evec1(5);
v1_6 = evec1(6);

v2_1 = evec2(1);
v2_2 = evec2(2);
v2_3 = evec2(3);
v2_4 = evec2(4);
v2_5 = evec2(5);
v2_6 = evec2(6);


v3_1 = evec3(1);
v3_2 = evec3(2);
v3_3 = evec3(3);
v3_4 = evec3(4);
v3_5 = evec3(5);
v3_6 = evec3(6);


v4_1 = evec4(1);
v4_2 = evec4(2);
v4_3 = evec4(3);
v4_4 = evec4(4);
v4_5 = evec4(5);
v4_6 = evec4(6);


v5_1 = evec5(1);
v5_2 = evec5(2);
v5_3 = evec5(3);
v5_4 = evec5(4);
v5_5 = evec5(5);
v5_6 = evec5(6);


v6_1 = evec6(1);
v6_2 = evec6(2);
v6_3 = evec6(3);
v6_4 = evec6(4);
v6_5 = evec6(5);
v6_6 = evec6(6);


%Load the algebraic solution to the expansion coefficients
load('wCoef_6state_condensed_solution.mat')


%Substitute the symbolic eigenvector components for their number in double
tic
c1_1 = double(subs(c1_1));
c1_2 = double(subs(c1_2));
c1_3 = double(subs(c1_3));
c1_4 = double(subs(c1_4));
c1_5 = double(subs(c1_5));
c1_6 = double(subs(c1_6));

c2_1 = double(subs(c2_1));
c2_2 = double(subs(c2_2));
c2_3 = double(subs(c2_3));
c2_4 = double(subs(c2_4));
c2_5 = double(subs(c2_5));
c2_6 = double(subs(c2_6));

c3_1 = double(subs(c3_1));
c3_2 = double(subs(c3_2));
c3_3 = double(subs(c3_3));
c3_4 = double(subs(c3_4));
c3_5 = double(subs(c3_5));
c3_6 = double(subs(c3_6));

c4_1 = double(subs(c4_1));
c4_2 = double(subs(c4_2));
c4_3 = double(subs(c4_3));
c4_4 = double(subs(c4_4));
c4_5 = double(subs(c4_5));
c4_6 = double(subs(c4_6));

c5_1 = double(subs(c5_1));
c5_2 = double(subs(c5_2));
c5_3 = double(subs(c5_3));
c5_4 = double(subs(c5_4));
c5_5 = double(subs(c5_5));
c5_6 = double(subs(c5_6));

c6_1 = double(subs(c6_1));
c6_2 = double(subs(c6_2));
c6_3 = double(subs(c6_3));
c6_4 = double(subs(c6_4));
c6_5 = double(subs(c6_5));
c6_6 = double(subs(c6_6));

elapsedTime = toc;
if clockMode == 1
    task_str = 'calculate cj_i by substituting in the eigenvector components';
    disp(['Took ' num2str(elapsedTime) ' seconds to ' task_str]);
end
%Takes about 3 seconds to sub in all the values


%-------------------------------------------------------------------------
% Calculate the conditional probabilities using matrices
%-------------------------------------------------------------------------

%Make a matrix C of expansion coefficient components
C =[c1_1, c1_2, c1_3, c1_4, c1_5, c1_6;...
   c2_1, c2_2, c2_3, c2_4, c2_5, c2_6;...
   c3_1, c3_2, c3_3, c3_4, c3_5, c3_6;...
   c4_1, c4_2, c4_3, c4_4, c4_5, c4_6;...
   c5_1, c5_2, c5_3, c5_4, c5_5, c5_6;...
   c6_1, c6_2, c6_3, c6_4, c6_5, c6_6];
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