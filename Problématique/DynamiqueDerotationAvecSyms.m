load('App6Problematique_Identification')
run constantes
run App6Problematique_Identification
syms Theta Thetadot

Delta_cmd_t1 = ((S_c*d_c*(C_Malpha*(Theta-Gamma_ref)+d_c./(2*V_RAA)*C_Mq*Thetadot))./C_Mdelta);
Delta_cmd_t2 = (400./((1/2)*Rho_0*exp(-H./H_s).*V_RAA.^2)*C_Mdelta);
Delta_cmd_t3 = 28./((S_c*d_c*(C_Malpha*(Theta-Gamma_ref)+d_c./(2*V_RAA)*C_Mq*Thetadot)).*((1/2)*Rho_0*exp(-H./H_s).*V_RAA.^2));
Delta_cmd    = Delta_cmd_t1+Delta_cmd_t2+Delta_cmd_t3;

