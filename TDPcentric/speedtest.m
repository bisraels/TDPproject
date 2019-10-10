syms a b c d x
 
r=solve(a*x^3+b*x^2+c*x+d);
 
R1 = inline(vectorize(char(r(1))),'a','b','c','d');
R2 = inline(vectorize(char(r(2))),'a','b','c','d');
R3 = inline(vectorize(char(r(3))),'a','b','c','d'); 
 
% enter some values to play with.
% note that there are 5 values for each coefficient.
atemp = rand(5,1)*3;
btemp = rand(5,1)*3;
ctemp = rand(5,1)*3;
dtemp = rand(5,1)*3;
 
% each roots variable will have 5 elements,
% which corresponds to each element in the coefficient
% vectors.
roots1 = R1(atemp,btemp,ctemp,dtemp);
roots2 = R2(atemp,btemp,ctemp,dtemp);
roots3 = R3(atemp,btemp,ctemp,dtemp);