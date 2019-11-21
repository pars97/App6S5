close all
clear all
clc

load('Accelero_Data_from_Moscow');
%% Définition paramètres et données connues
Temps        = t;           % Temps des données des russes en s
Acc_mes      = acc_mes;     % Données des russes en m/s^2
M_c          = 50;          % Masse de la capsule en kg
J_c          = 1.5;         % Inertie de la capsule en kg-m^2
R_mars       = 3397*10^3;   % Rayon de mars en m
Mu_mars      = 42830*10^9;  % Parametre gravitationelle de mars en m^3/s^2
S_c          = 0.8;         % Surface aero. de la capsule en m^2
d_c          = 0.05;        % Dimension de la capsule en m
C_do         = 1.2;         % Coefficient de trainée
C_Lalpha     = 0.8;         % Coefficient de portance
C_Malpha     = -0.07;       % Coefficient de couple
C_Mq         = -0.05;       % Coefficient d'amortissement
C_Mdelta     = 0.1;         % Coefficient de volet aero.
V_IniR       = 6100;        % Vitesse initiale expérience des russes en m/s
Gamma_IniR   = -90;         % Angle gamma initial expérience des russes en °
H_IniR       = 120000;      % Hauteur initiale expérience des russes en m
Theta_IniR   = -90;         % Angle theta initial expérience des russes en °
Br_Mes_AccR  = 0.03;        % Bruit sur les mesur de l'accéléromètre des russes
P_BalR       = 0.0192;      % Paramètre balistique des russe en m^2/kg 
V_IniM       = 6100;        % Vitesse initiale du mandat en m/s
Gamma_IniM   = -20.5;       % Angle gamma initial du mandat en °
H_IniM       = 120000;      % Hauteur initiale du mandat en m
Theta_IniM   = -80;         % Angle theta initial du mandat en °
V_FinM1      = 250;         % Vitesse maximale optimale finale du mandat en m/s
V_FinM2      = 300;         % Vitesse maximale finale du mandat en m/s
H_FinM       = 10000;       % Hauteur du déploiment des parachutes (h finale) en m
T_Lim        = 45;          % Temps limite pour une force sur la capsule de plus de 2000N en s
Theta_cmdLim = 60;          % Angle maximale en absolue de la capsule en °

%% Détermination des paramètre Hs et Po
H_V_Trap = Temps(2)-Temps(1);

% Détermination de la vitesse selon l'accélération par méthode des trapèze
DpAccMes_a = (Acc_mes(2)-Acc_mes(1))/H_V_Trap;
DpAccMes_b = (Acc_mes(end)- Acc_mes(end-1))/H_V_Trap;

V_Mes_Trap   = H_V_Trap/2*(Acc_mes(1)+Acc_mes(end)+2*sum(Acc_mes(2:end-1)))+V_IniM;
E_V_Mes_Trap = -H_V_Trap^2/12*(DpAccMes_b-DpAccMes_a);

%intégrale point par point
V_Mes_Trap_Vec(1)=V_IniM;
for n=2:1:length(Temps)
    V_Mes_Trap_Vec(n) = H_V_Trap/2*(Acc_mes(1)+Acc_mes(n)+2*sum(Acc_mes(2:n-1)))+V_IniM;
end


% Détermination de la vitesse selon l'accélération par méthode de Simpson
DtAccMes_a = (Acc_mes(4)-3*Acc_mes(3)+3*Acc_mes(2)-Acc_mes(1))/H_V_Trap^3;
DtAccMes_b = (Acc_mes(end)-3*Acc_mes(end-1)+3*Acc_mes(end-2)-Acc_mes(end-3))/H_V_Trap^3;

V_Mes_Simp = H_V_Trap/3*(Acc_mes(1)+Acc_mes(end)+4*sum(Acc_mes(2:2:end-1))+2*sum(Acc_mes(3:2:end-1)))+V_IniM;
E_V_Mes_Simp = -H_V_Trap^4/180*(Acc_mes(end)-Acc_mes(1))*(DtAccMes_b-DtAccMes_a);

V_Mes_Simp_Aff(1)=V_IniM;
for n=3:2:length(Temps)
    V_Mes_Simp_Aff(((n-1)/2)+1) = H_V_Trap/3*(Acc_mes(1)+Acc_mes(n)+4*sum(Acc_mes(2:2:n-1))+2*sum(Acc_mes(3:2:n-1)))+V_IniM;
end
figure()
plot (Temps(1:2:end),V_Mes_Simp_Aff)
hold on
plot(Temps,V_Mes_Trap_Vec);
title('Vitesse en fonction du temps selon différente approximation')
legend('Simpson','Trapèze')
% La méthode des trapèze est sélectionnée car sont erreure est plus petite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Détermination de la hauteur selon la vitesse
DpVmes_a = (V_Mes_Trap_Vec(2)-V_Mes_Trap_Vec(1))/H_V_Trap;
DpVmes_b = (V_Mes_Trap_Vec(end)- V_Mes_Trap_Vec(end-1))/H_V_Trap;

H_Mes_Trap   = H_IniM - H_V_Trap/2*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(end)+2*sum(V_Mes_Trap_Vec(2:end-1)));
E_H_Mes_Trap = -H_V_Trap^2/12*(DpVmes_b-DpVmes_a)

%intégrale point par point
H_Mes_Trap_Aff(1)=H_IniM;
for n=2:1:length(Temps)
    H_Mes_Trap_Aff(n) = H_IniM - H_V_Trap/2*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(n)+2*sum(V_Mes_Trap_Vec(2:n-1)));
end

% Détermination de la vitesse selon l'accélération par méthode de Simpson
DtHmes_a = (V_Mes_Trap_Vec(4)-3*V_Mes_Trap_Vec(3)+3*V_Mes_Trap_Vec(2)-V_Mes_Trap_Vec(1))/H_V_Trap^3;
DtHmes_b = (V_Mes_Trap_Vec(end)-3*V_Mes_Trap_Vec(end-1)+3*V_Mes_Trap_Vec(end-2)-V_Mes_Trap_Vec(end-3))/H_V_Trap^3;

H_Mes_Simp = H_IniM - H_V_Trap/3*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(end)+4*sum(V_Mes_Trap_Vec(2:2:end-1))+2*sum(V_Mes_Trap_Vec(3:2:end-1)));
E_H_Mes_Simp = -H_V_Trap^4/180*(V_Mes_Trap_Vec(end)-V_Mes_Trap_Vec(1))*(DtHmes_b-DtHmes_a)

H_Mes_Simp_Aff(1)=H_IniM;
for n=3:2:length(Temps)
    H_Mes_Simp_Aff(((n-1)/2)+1) = H_IniM - H_V_Trap/3*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(n)+4*sum(V_Mes_Trap_Vec(2:2:n-1))+2*sum(V_Mes_Trap_Vec(3:2:n-1)));
end
figure()
plot (Temps(1:2:end),H_Mes_Simp_Aff)
hold on
plot(Temps,H_Mes_Trap_Aff);
% title('Hauteur en fonction du temps selon différente approximation')
legend('Simpson','Trapèze')