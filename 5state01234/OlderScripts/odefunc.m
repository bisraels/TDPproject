% function for odefunc

function dydt = odefunc(t,y)
dydt = zeros(5,1);
% dydt(1) = y(1)+2*y(2);
% dydt(2) = 3*y(1)+2*y(2);
 M = [1, 2,3,4,5; 3, 2,1,4,5;1,5,2,4,3;5,4,3,2,1;5,1,4,2,3];
%M = [1,2;4,5];
dydt = diff(y) == M * y;
% dydt = M * y;

end

% function dPdt = odefunc(t,P,K,k01,k02,k03,k04,k10,k12,k13,k14,k20,k21,k23,k24,k30,k31,k32,k34,k40,k41,k42,k43)
% dPdt = zeros(5,1);
% % dydt(1) = y(1)+2*y(2);
% % dydt(2) = 3*y(1)+2*y(2);
% % P(t) = [P0(t); P1(t); P2(t); P3(t); P4(t)];
% 
% % k32 = (k23 * k34 * k42)/(k24 * k43);
% % 
% % P = [P0(t); P1(t); P2(t); P3(t); P4(t)];
% 
% k01 = 0.001;
% k02 = 2;
% k03 = 0;
% k04 = 0;
% k10 = 0.000034;
% k12 = 0;
% k13 = 0;
% k14 = 0;
% k20 = 0.153;
% k21 = 0;
% k23 = 700;
% k24 = 34;
% k30 = 0; 
% k31 = 0;
% % k32 = 87;
% k34 = 45;
% k40 = 0;
% k41 = 0;
% k42 = 4;
% k43 = 7;
% k32 = (k23 * k34 * k42)/(k24 * k43);
% 
% K = [-(k01+k02+k03+k04), k10,             k20,              k30,              k40;...
%         k01,      -(k10+k12+k13+k14),     k21,              k31,              k41;...
%         k02,             k12,       -(k20+k21+k23+k24),     k32,              k42;...
%         k03,             k13,             k23,      -(k30+k31+k32+k34),       k43;...
%         k04,             k14,             k24,              k34,        -(k40+k41+k42+k43)];
% 
% P = [P0(t); P1(t); P2(t); P3(t); P4(t)];
% 
% 
% dPdt = diff(P)==K * P;
% %dPdt = K * P;
% 
% end
