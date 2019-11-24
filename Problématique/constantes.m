%% constantes
addpath ../
load('Accelero_Data_from_Moscow');
%% Définition paramètres et données connues
Temps        = t;                   % Temps des données des russes en s
Acc_mes      = -acc_mes;            % Données des russes en m/s^2
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
V_IniR       = 6100;                % Vitesse initiale expérience des russes en m/s
Gamma_IniR   = deg2rad(-90);        % Angle gamma initial expérience des russes en rad
H_IniR       = 120000;              % Hauteur initiale expérience des russes en m
Theta_IniR   = deg2rad(-90);        % Angle theta initial expérience des russes en °
Br_Mes_AccR  = 0.03;                % Bruit sur les mesur de l'accéléromètre des russes
P_BalR       = 0.0192;              % Paramètre balistique des russe en m^2/kg 
V_IniM       = 6100;                % Vitesse initiale du mandat en m/s
Gamma_IniM   = deg2rad(-20.5);      % Angle gamma initial du mandat en °
H_IniM       = 120000;              % Hauteur initiale du mandat en m
Theta_IniM   = deg2rad(-80);        % Angle theta initial du mandat en °
V_FinM       = [250 300]';           % Described below
V_FinM1      = 250;                 % Vitesse maximale optimale finale du mandat en m/s
V_FinM2      = 300;                 % Vitesse maximale finale du mandat en m/s
H_FinM       = 10000;               % Hauteur du déploiment des parachutes (h finale) en m
T_Lim        = 45;                  % Temps limite pour une force sur la capsule de plus de 2000N en s
Theta_cmdLim = deg2rad(60);         % Angle maximale en absolue de la capsule en °
B            = S_c*C_do/M_c;        % Paramètre Balistique constant (m^2/kg)
disp('Constantes loaded')