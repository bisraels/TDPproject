% FourStateODE
%
% y(1) -> 0 bound, y(2) -> 1 bound and no transition to 2 bound, y(3) -> 1
% bound and transition possiblity to 2 bound, y(4) -> 2 bound.

function dydt = FourStateODE(t,y,opt,k01,k10,k02,k20,k03,k30,k12,k21,k13,k31,k23,k32)
dydt = zeros(size(y));

dydt(1) = -(k01 + k02 + k03)*y(1) + k10*y(2) + k20*y(3) + k30*y(4);
dydt(2) = k01*y(1) - (k10 + k12 + k13)*y(2) + k21*y(3) + k31*y(4);
dydt(3) = k02*y(1) + k12*y(2) - (k20 + k21 + k23)*y(3) + k32*y(4);
dydt(4) = k03*y(1) + k13*y(2) + k23*y(3) - (k30 + k31 + k32)*y(4);

end