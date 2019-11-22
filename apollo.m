function f = apollo(t,VE)


u = 1/82.45;
u_star = 1-u;

r1 = sqrt((VE(1)+u)^2+VE(2)^2);
r2 = sqrt((VE(1)-u_star)^2+VE(2)^2);

f(3) = 2*VE(4)+VE(1)-u_star*(VE(1)+u)/r1^3-u*(VE(1)-u_star)/r2^3;
f(4) = -2*VE(3)+VE(2)-u_star*VE(2)/r1^3-u*VE(2)/r2^3;
f(1) = VE(3);
f(2) = VE(4);

  f = f(:);