function f = ass1(t, VE)

Rho_0 = 0.021571989401399;
H_s = 1.100353042442160e+04;

ctrl =1;
tau = 0.2;
zeta = 0.7;
wn = 20;

    %% Définition paramètres et données connues
M_c          = 50;                  % Masse de la capsule en kg
J_c          = 1.5;                 % Inertie de la capsule en kg-m^2
R_mars       = 3397*10^3;           % Rayon de mars en m
Mu_mars      = 42830*10^9;          % Parametre gravitationelle de mars en m^3/s^2
S_c          = 0.8;                 % Surface aero. de la capsule en m^2
d_c          = 0.05;                % Dimension de la capsule en m
C_do         = 1.2;                 % Coefficient de trainée
C_Lalpha     = 0.8;                 % Coefficient de portance
C_Malpha     = -0.07;               % Coefficient de couple
C_Mq         = -0.05;               % Coefficient d'amortissement
C_Mdelta     = 0.1;                 % Coefficient de volet aero.
V_FinM1      = 250;                 % Vitesse maximale optimale finale du mandat en m/s
V_FinM2      = 300;                 % Vitesse maximale finale du mandat en m/s
H_FinM       = 10000;               % Hauteur du déploiment des parachutes (h finale) en m
Theta_cmdLim = deg2rad(60);         % Angle maximale en absolue de la capsule en °
B            = S_c*C_do/M_c;        % Paramètre Balistique constant (m^2/kg)
   
%%

    Rho = Rho_0*exp(-VE(3)/H_s);
    
    r = R_mars + VE(3);
    
    r_fin = R_mars + H_FinM;
    
    g = Mu_mars / r.^2;
    
    alpha = VE(5) - VE(2);
    
    Rho_fin = Rho_0*exp(-H_FinM/H_s);
    
    Pdyn = (1/2) * Rho * (VE(1).^2);
    
    B = C_do * S_c / M_c;
    
    DVA = V_FinM2 - sqrt(VE(1)^2+((2*Mu_mars)*((1/r_fin)-(1/r))));

    gamma_ref = asin((1/2)*B*H_s*((Rho_fin - Rho) ./ log(1 + DVA/VE(1))));
        
    g_gamma = (Pdyn * S_c * C_Lalpha) / (VE(1) * M_c); 
         
    theta_eq = VE(2) - (cos(VE(2)) * M_c) / (Pdyn*S_c*C_Lalpha) * (VE(1)^2/r - Mu_mars / r^2);
        
    theta_cmd = theta_eq + (1/(tau*g_gamma)) * (gamma_ref - VE(2));
  
    if (theta_cmd > deg2rad(60))
            theta_cmd = deg2rad(60);
        
    end
    
    if (theta_cmd < -deg2rad(60))
            theta_cmd = -deg2rad(60);
        
    end
    
    
    Kp = wn^2;
        
    Kd = 2 * zeta * wn;
        
    k = (C_Mdelta * Pdyn * S_c * d_c) / J_c;
        
    delta_eq = (-C_Malpha/C_Mdelta) * (VE(5) - VE(2)) - (C_Mq * d_c * VE(6)) / (2 * VE(1) * C_Mdelta);
        
    delta_cmd = delta_eq + Kp/k * (theta_cmd - VE(5)) + Kd/k * (-VE(6));
 
    Daero = Pdyn * S_c * C_do;
    
    Laero = Pdyn * S_c * C_Lalpha * alpha;
    
    Maero = Pdyn * S_c * d_c * (C_Malpha * alpha + (d_c / (2*VE(1))) * C_Mq*VE(6) + C_Mdelta*delta_cmd); 
    
    f(1) = - Daero / M_c - g*sin(VE(2));
    f(2) = (1/VE(1)) * (Laero/M_c + ((VE(1).^2/r) - g)*cos(VE(2)));
    f(3) = VE(1)*sin(VE(2));
    f(4) = (VE(1)/r) * cos(VE(2));
    f(5) = VE(6);
    f(6) = (1/J_c) * Maero;
    
    
    f = f(:);
end


