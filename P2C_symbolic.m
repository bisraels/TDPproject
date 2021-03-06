% function P2C2()

%P will be the NxNxt matrix. Default is a 2x2 matrix
P = k2P();

[N,~,timesteps] = size(P);

% Make a symbolic column vector of the discrete probability
Peq = sym('P%d_eq',[N 1],'positive')  ; %Peq(i) = Pi_eq

%Make a symbolic row vector of the FRET Values
A = sym('A%d',[1 N],'positive');         %A(i) = Ai

cP = sym('P',[N N]);        %cP(m,n) = Pm_n

% Define symbolic matrix of c's for our constants
c = sym('c', [N N]);                    %c(m,n) = cm_n

% Define symbolic matrix of the eigenvector components.
% Transpose the matrix for proper matrix mulitplication
v = sym('v',[N N]).';                   %v(m,n) = vn_m

% Make a symbolic column vector of the eigenvalues and time
lam = exp(sym('eval',[N 1])*'time');    %lam(i) = exp(evali*time)

% % %make a symbolic figure
% % p = v*(c.*lam);

% p(:)

% C2_sim = P3_eq*A3*(A3*p3_3 + A2*p2_3 + A1*p1_3) +...
%     P2_eq*A2*(A3*p3_2 + A2*p2_2 + A1*p1_2) +...
%     P1_eq*A1*(A3*p3_1 + A2*p2_1 + A1*p1_1);

%--------------------------------------------------------------------------
% Calculate 2 point TCF (Sums using loops)
%--------------------------------------------------------------------------

if exist('C2','var')
    clear('C2');
end
syms C2;
% C2 = [];
for i = 1:numel(A)
    for j = 1:numel(A)
        if i == 1 && j == 1
%                         C2 = A(j)*p(j,i)*A(i)*Peq(i);   %Long form with eigenvectors, eigenvalues, and expansion coeficcients
            C2 = A(j)*cP(j,i)*A(i)*Peq(i); %short form with just P3_3
        else
%                         C2temp = A(j)*p(j,i)*A(i)*Peq(i);
            C2temp = A(j)*cP(j,i)*A(i)*Peq(i);
            C2 =  C2 + C2temp;
        end
    end
end

%For each row/column, 1:N, set the conditional probabilites equal.
 if N >= 2
    P1_1 = reshape(P(1,1,:),[1 timesteps]);
    P1_2 = reshape(P(1,2,:),[1 timesteps]);
    P2_1 = reshape(P(2,1,:),[1 timesteps]);
    A1 = 0.1;
    A2 = 0.2;
 end
 
 
%--------------------------------------------------------------------------
% Calculate 2 point TCF (Matrix)
%--------------------------------------------------------------------------
% 
% Amat = A.*eye(numel(A));
% % C2test = Amat*p*Amat*Peq;
% C2test = sum(Amat*cP*Amat*Peq);

% return
%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------

syms C4
for i = 1:numel(A)
    for j = 1:numel(A)
        for k = 1:numel(A)
            for l = 1:numel(A)
                if i == 1 && j == 1 && k == 1 && l == 1
                    %             C4 = A(l)*p(l,k)*A(k)*p(k,j)*A(j)*p(j,i)*A(i)*Peq(i);  %Long form with eigenvectors, eigenvalues, and expansion coeficcients
                    C4 = A(l)*cP(l,k)*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i); %short form with just P3_3
                else
                    %             C4temp = A(l)*p(l,k)*A(k)*p(k,j)*A(j)*p(j,i)*A(i)*Peq(i);
                    C4temp = A(l)*cP(l,k)*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i);
                    C4 =  C4 + C4temp;
                end
            end
        end
    end
end