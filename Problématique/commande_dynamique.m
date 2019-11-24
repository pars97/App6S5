%% Commandes de la dynamique
clear all
close all
clc
run constantes
load ('App6Problematique_Identification.mat')
rfin = (R_mars+H_FinM);
Rho_fin = Rho_0*exp(-H_FinM/H_s);
Rho_Ini = Rho_0*exp(-H_IniM/H_s);


Rho = Rho_0*exp(-H/H_s);
r = (R_mars+H);

% V(1,1) = V_IniM*exp(1/2*B*H_s*(Rho(2)-Rho(1))/sin(Gamma_IniM)); 
% V(2,1) = V_IniM*exp(1/2*B*H_s*(Rho(2)-Rho(1))/sin(Gamma_IniM)); 
% Delta_V_aero(:,1) = V_FinM-sqrt(V(1)^2+2*Mu_mars*(1/rfin-1/r(1)));
% Gamma_ref(:,1) = asin(1/2*B*H_s*(Rho_fin-Rho(1))./log(1+Delta_V_aero(:,1)./V(:,1)));
%%


V(1) = V_IniM*exp(1/2*B*H_s*(Rho(1)-Rho_Ini)/sin(Gamma_IniM)); 
Delta_V_aero(1) = V_FinM(1)-sqrt(V(1)^2+2*Mu_mars*(1/rfin-1/r(1)));
Gamma_ref(1) = asin(1/2*B*H_s*(Rho_fin-Rho(1))./log(1+Delta_V_aero(1)/V(1)));


for i = 1:1:length(H)-1

V(i+1) = V_IniM*exp(1/2*B*H_s.*(Rho(i)-Rho(1))./sin(Gamma_ref(i)));

Delta_V_aero(i+1) = V_FinM(1)-sqrt(V(i)^2+2*Mu_mars.*(1/rfin-1/r(i)));
Gamma_ref(i+1) = asin(1/2*B*H_s*(Rho_fin-Rho(i+1))/(log(1+Delta_V_aero(i+1)/V(i+1))));


end



% 


