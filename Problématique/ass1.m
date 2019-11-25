function f = ass1(t,VE)
run constantes
load ('App6Problematique_Identification.mat')


if (VE(5)>deg2rad(60))
    VE(5) = deg2rad(60);
end
if (VE(5)<deg2rad(-60))
    VE(5) = deg2rad(-60);
end

Rho = Rho_0*exp(-VE(3)/H_s);
Rho_dyn = 1/2*Rho*VE(1)^2;
alpha = VE(5)-VE(2);
Daero = Rho_dyn*S_c*C_do;
Laero = Rho_dyn*S_c*C_Lalpha*alpha;
r = R_mars + VE(3);
gr = Mu_mars/r^2;

DVA = V_FinM2 - sqrt(VE(1)^2+2*Mu_mars*(1/(R_mars+H_FinM)-1/r));
Gamma_ref = asin(1/2*B*H_s*(Rho_0*exp(-H_FinM/H_s)-Rho)/log(1+DVA/VE(1)));


f_gamma = -Rho_dyn*S_c*C_Lalpha*VE(2)/(M_c*VE(1))+cos(VE(2))/VE(1)*(VE(1)^2/r-gr);
g_gamma = Rho_dyn*S_c*C_Lalpha*VE(5)/(M_c*VE(1));

Kp_gamma =5;


Theta_cmd = -f_gamma/g_gamma+Kp_gamma/g_gamma*(Gamma_ref-VE(2));


Delta_cmd_t1 = -(((S_c*d_c*(C_Malpha*alpha)+d_c./(2*VE(1))*C_Mq*VE(6)))/C_Mdelta);
Delta_cmd_t2 = (400/((1/2)*Rho_0*exp(-VE(3)/H_s)*VE(1)^2)*C_Mdelta)*(Theta_cmd-VE(5));
Delta_cmd_t3 = (28/(((S_c*d_c*(C_Malpha*alpha)+d_c/(2*VE(1))*C_Mq*VE(6)))*((1/2)*Rho_0*exp(-VE(3)/H_s)*VE(1)^2)))*(0-VE(6));
Delta_cmd    = Delta_cmd_t1+Delta_cmd_t2+Delta_cmd_t3;
% Delta_cmd =0;


Maero = Rho_dyn*S_c*d_c*(C_Malpha*alpha+d_c/(2*VE(1))*C_Mq*VE(6)+C_Mdelta*Delta_cmd);





f(1) = -Daero/M_c-gr*sin(VE(2));
f(2) = f_gamma+g_gamma*(-f_gamma/g_gamma + Kp_gamma/g_gamma*(Gamma_ref-VE(2)));
f(3) = VE(1)*sin(VE(2));
f(4) = VE(1)/r*cos(VE(2));
f(5) = VE(6);
f(6) = 1/J_c*Maero;
  f = f(:);


