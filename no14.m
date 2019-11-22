function f = no14(t,VE)

xcmd = 10;
vxcmd = 0;

Wn = 4;
zeta = 0.5;

Kp = Wn^2;
Kd = 2*zeta*Wn;

u = 0.6*VE(2)+3*VE(1)+VE(1)^2 + Kp*(xcmd-VE(1)) + Kd*(vxcmd-VE(2));
f(1) = VE(2);
f(2) = -0.6*VE(2)-3*VE(1)-VE(1)^2+u;

  f = f(:);