% function P2C2()

% switch nargin
%     case 0
     [P] = k2P();
% end
% 
[N,~,length_time] = size(P);


% Make a symbolic column vector of the discrete probability
% Peq = sym('P%d_eq',[N 1],'positive')  ; %Peq(i) = Pi_eq
Peq = P(:,1,end);   %The first element of each row is Pj,1(inf) (eq pop)

%Make a symbolic row vector of the FRET Values
% A = sym('A%d',[1 N],'positive');         %A(i) = Ai
A = linspace(0,1,N);

% Make a symbolic matrix of the conditional probability
% cP = sym('P',[N N]);        %cP(m,n) = Pm_n


%--------------------------------------------------------------------------
% Calculate 2 point TCF (with loop)
%--------------------------------------------------------------------------
time = time_sim;
Npts = numel(time);
C2_sim = zeros(size(time));
%Subtract the mean FRET 
Amean = sum(Peq.*A');
A = A - Amean;
tic
for i = 1:numel(A)
    for j = 1:numel(A)
        C2_sim_temp = A(j) * reshape(P(j,i,:),[1,Npts]) * A(i) * Peq(i);
        C2_sim = C2_sim + C2_sim_temp;
    end
end
elapsedTime = toc;
disp(['Time to calculate C2 for N = ' num2str(N) ' = ' num2str(elapsedTime) ' seconds']);



%--------------------------------------------------------------------------
% Calculate 2 point TCF (with matrix)
%--------------------------------------------------------------------------
%%
A = linspace(0,1,N);
Amat = A.*eye(N);
Amat3D = ones(N,N,length(time));
Peq2D = ones(N,length(time));
for i = 1:length(time)
    Amat3D(:,:,i) = Amat;
    Peq2D(:,i) = Peq;
end
% C2test = sum(Amat*P*Amat*Peq);
C2test = sum(Amat3D*P*Amat3D*Peq2D);



%--------------------------------------------------------------------------
% Calculate 4 point TCF
%--------------------------------------------------------------------------
% 
% syms C4
% for i = 1:numel(A)
%     for j = 1:numel(A)
%         for k = 1:numel(A)
%             for l = 1:numel(A)
%                 if i == 1 && j == 1 && k == 1 && l == 1
%                     %             C4 = A(l)*p(l,k)*A(k)*p(k,j)*A(j)*p(j,i)*A(i)*Peq(i);  %Long form with eigenvectors, eigenvalues, and expansion coeficcients
%                     C4 = A(l)*cP(l,k).'*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i); %short form with just P3_3
%                 else
%                     %             C4temp = A(l)*p(l,k)*A(k)*p(k,j)*A(j)*p(j,i)*A(i)*Peq(i);
%                     C4temp = A(l)*cP(l,k).'*A(k)*cP(k,j)*A(j)*cP(j,i)*A(i)*Peq(i);
%                     C4 =  C4 + C4temp;
%                 end
%             end
%         end
%     end
% end