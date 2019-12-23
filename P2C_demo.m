% function P2C2()

% [N,~,length_time] = size(P);
N = 2;

% Make a symbolic column vector of the discrete probability
Peq = sym('P%d_eq',[N 1],'positive')  ; %Peq(i) = Pi_eq

%Make a symbolic row vector of the FRET Values
A = sym('A%d',[1 N],'positive');         %A(i) = Ai

% p(:)

%--------------------------------------------------------------------------
% Calculate 2 point TCF
%--------------------------------------------------------------------------

if exist('C2','var')
    clear('C2');
end
syms C2;
% Make a symbolic matrix of the conditional probability
cP = sym('P',[N N]);        %cP(m,n) = Pm_n

% C2 = [];
for i = 1:numel(A)
    for j = 1:numel(A)
        if i == 1 && j == 1
            %             C2 = A(j)*p(j,i)*A(i)*Peq(i);   %Long form with eigenvectors, eigenvalues, and expansion coeficcients
            C2 = A(j)*cP(j,i)*A(i)*Peq(i); %short form with just P3_3

        else
            %             C2temp = A(j)*p(j,i)*A(i)*Peq(i);
            C2temp = A(j)*cP(j,i)*A(i)*Peq(i);
            C2 =  C2 + C2temp;
        end
    end
end

%--------------------------------------------------------------------------
% Calculate 2 point TCF (Matrix)
%--------------------------------------------------------------------------

% Amat = A.*eye(N);
% C2test = sum(Amat*cP*Amat*Peq);
% time = time_sim;
% Amat3D = Amat;
% Peq2D = Peq;
% for i = 1:length(time)-1
%     Amat3D = cat(3,Amat3D,Amat);
%     Peq2D(:,i) = cat(Peq2D,Peq);
% end

% C2 = Amat3D
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
                    C4 = A(l)*cP(l,k).'*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i); %short form with just P3_3
                else
                    %             C4temp = A(l)*p(l,k)*A(k)*p(k,j)*A(j)*p(j,i)*A(i)*Peq(i);
                    C4temp = A(l)*cP(l,k).'*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i);
                    C4 =  C4 + C4temp;
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% Solve for conditional probabilites using the eigen basis
%--------------------------------------------------------------------------


% Define symbolic matrix of c's for our constants
c = sym('c', [N N]);                    %c(m,n) = cm_n

% Define symbolic matrix of the eigenvector components.
% Transpose the matrix for proper matrix mulitplication
v = sym('v',[N N]).';                   %v(m,n) = vn_m

% Make a symbolic column vector of the eigenvalues and time
lam = exp(sym('eval',[N 1])*'time');    %lam(i) = exp(evali*time)

%make a symbolic figure
p = v*(c.*lam);

k = sym('k%d%d',[N,N]).';
for i = 1:N
    k(i,i) = 0;
    k(i,i) = -sum(k(:,i));
end


