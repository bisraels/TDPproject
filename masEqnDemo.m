% symbolicDemo

k = sym('k',[N N]).'
for i = 1:N
    k(i,i) = 0;
    k(i,i) = -sum(k(:,i));
end

Pdot = k*P