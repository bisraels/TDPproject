
%Number of states
N = 6

% Define symbolic matrix of c's for our constants
c = sym('c', [N,N])

% Define symbolic matrix of the eigenvector components
v = sym('v',[N,N]).' 

% Make a symbolic vector of the eigenvalues and time
lam = exp(sym('lam',[N,1])*'time')

% time = 1:10;
p = v*(c.*lam)


%Simulate a K matrix for a 6-state system
K = [-25,21,31,0,0,0;...
    12,-93,2852/91,42,52,0;...
    13,23,-(31+2852/91),0,0,0;...
    0,24,0,-42,0,0;...
    0,25,0,0,-108,65;...
    0,0,0,0,56,-65];

sum(K);

[Vec, lambda] = eig(K);

[lambdaSort, index] = sort(diag(lambda),'descend');   % sort just in case, it the program seems to output, closest to 0 -> most negative
lambdaSorted = lambda(index,index);
VecSorted = Vec(:,index);


lam1 = double(lambdaSorted(1,1));
lam2 = double(lambdaSorted(2,2));
lam3 = double(lambdaSorted(3,3));
lam4 = double(lambdaSorted(4,4));
lam5 = double(lambdaSorted(5,5));
lam6 = double(lambdaSorted(6,6));

evec1 = VecSorted(:,1);
evec2 = VecSorted(:,2);
evec3 = VecSorted(:,3);
evec4 = VecSorted(:,4);
evec5 = VecSorted(:,5);
evec6 = VecSorted(:,6);


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
tic
%Substitute the symbolic eigenvector components for their number in double
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

toc
C =[c1_1, c1_2, c1_3, c1_4, c1_5, c1_6;...
   c2_1, c2_2, c2_3, c2_4, c2_5, c2_6;...
   c3_1, c3_2, c3_3, c3_4, c3_5, c3_6;...
   c4_1, c4_2, c4_3, c4_4, c4_5, c4_6;...
   c5_1, c5_2, c5_3, c5_4, c5_5, c5_6;...
   c6_1, c6_2, c6_3, c6_4, c6_5, c6_6];
disp(C)

V = VecSorted.';

N = 6;
Npts = 150;
time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;

lam = exp(sym('lam',[1,N])*'time')
Lam1 = double(subs(lam(1)));
Lam2 = double(subs(lam(2)));
Lam3 = double(subs(lam(3)));
Lam4 = double(subs(lam(4)));
Lam5 = double(subs(lam(5)));
Lam6 = double(subs(lam(6)));
LamVec = [Lam1;Lam2;Lam3;Lam4;Lam5;Lam6];

Npts = 150;
time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;


% time = 1:10;
% P = V*(C.*LamVec)
