close all
clear all
clc

run constantes
%% Détermination des paramètre Hs et Po
H_V_Trap = Temps(2)-Temps(1);

% Détermination de la vitesse selon l'accélération par méthode des trapèze
DpAccMes_a = (Acc_mes(2)-Acc_mes(1))/H_V_Trap;
DpAccMes_b = (Acc_mes(end)- Acc_mes(end-1))/H_V_Trap;


V_Mes_Trap   = H_V_Trap/2*(Acc_mes(1)+Acc_mes(end)+2*sum(Acc_mes(2:end-1)))+V_IniR;
E_V_Mes_Trap = -H_V_Trap^2/12*(DpAccMes_b-DpAccMes_a);

%intégrale point par point
V_Mes_Trap_Vec(1)=V_IniR;


for n=2:1:length(Temps)
    V_Mes_Trap_Vec(n) = H_V_Trap/2*(Acc_mes(1)+Acc_mes(n)+2*sum(Acc_mes(2:n-1)))+V_IniR;
end


% Détermination de la vitesse selon l'accélération par méthode de Simpson
DtAccMes_a = (Acc_mes(4)-3*Acc_mes(3)+3*Acc_mes(2)-Acc_mes(1))/H_V_Trap^3;
DtAccMes_b = (Acc_mes(end)-3*Acc_mes(end-1)+3*Acc_mes(end-2)-Acc_mes(end-3))/H_V_Trap^3;

V_Mes_Simp = H_V_Trap/3*(Acc_mes(1)+Acc_mes(end)+4*sum(Acc_mes(2:2:end-1))+2*sum(Acc_mes(3:2:end-1)))+V_IniR;
E_V_Mes_Simp = -H_V_Trap^4/180*(Acc_mes(end)-Acc_mes(1))*(DtAccMes_b-DtAccMes_a);

V_Mes_Simp_Vec(1)=V_IniR;



for n=3:2:length(Temps)
    V_Mes_Simp_Vec(((n-1)/2)+1) = H_V_Trap/3*(Acc_mes(1)+Acc_mes(n)+4*sum(Acc_mes(2:2:n-1))+2*sum(Acc_mes(3:2:n-1)))+V_IniR;
end



% figure()
% plot (Temps(1:2:end),V_Mes_Simp_Vec)
% hold on
% plot(Temps,V_Mes_Trap_Vec);
% title('Vitesse en fonction du temps selon différente approximation')
% legend('Simpson','Trapèze')
% La méthode des trapèze est sélectionnée car sont erreure est plus petite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Détermination de la hauteur selon la vitesse
DpVmes_a = (V_Mes_Trap_Vec(2)-V_Mes_Trap_Vec(1))/H_V_Trap;
DpVmes_b = (V_Mes_Trap_Vec(end)- V_Mes_Trap_Vec(end-1))/H_V_Trap;

H_Mes_Trap   = H_IniR + H_V_Trap/2*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(end)+2*sum(V_Mes_Trap_Vec(2:end-1)));
E_H_Mes_Trap = H_V_Trap^2/12*(DpVmes_b-DpVmes_a);

%intégrale point par point
H_Mes_Trap_Vec(1)=H_IniR;


for n=2:1:length(Temps)
    H_Mes_Trap_Vec(n) = H_IniR - H_V_Trap/2*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(n)+2*sum(V_Mes_Trap_Vec(2:n-1)));
end

% Détermination de la vitesse selon l'accélération par méthode de Simpson
DtHmes_a = (V_Mes_Trap_Vec(4)-3*V_Mes_Trap_Vec(3)+3*V_Mes_Trap_Vec(2)-V_Mes_Trap_Vec(1))/H_V_Trap^3;
DtHmes_b = (V_Mes_Trap_Vec(end)-3*V_Mes_Trap_Vec(end-1)+3*V_Mes_Trap_Vec(end-2)-V_Mes_Trap_Vec(end-3))/H_V_Trap^3;

H_Mes_Simp = H_IniR + H_V_Trap/3*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(end)+4*sum(V_Mes_Trap_Vec(2:2:end-1))+2*sum(V_Mes_Trap_Vec(3:2:end-1)));
E_H_Mes_Simp = -H_V_Trap^4/180*(V_Mes_Trap_Vec(end)-V_Mes_Trap_Vec(1))*(DtHmes_b-DtHmes_a);

H_Mes_Simp_Aff(1)=H_IniR;


for n=3:2:length(Temps)
    H_Mes_Simp_Aff(((n-1)/2)+1) = H_IniR - H_V_Trap/3*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(n)+4*sum(V_Mes_Trap_Vec(2:2:n-1))+2*sum(V_Mes_Trap_Vec(3:2:n-1)));
end


% figure()
% plot (Temps(1:2:end),H_Mes_Simp_Aff)
% hold on
% plot(Temps,H_Mes_Trap_Vec);
% % title('Hauteur en fonction du temps selon différente approximation')
% legend('Simpson','Trapèze')
% La méthode des trapèzes est sélectionnée car sont erreure est plus petite

%% Identification des paramètres et RMS/R2 linéaire

Rho_mes = (acc_mes'*2*M_c)./(V_Mes_Trap_Vec.^2*S_c*C_do);
Y = log(Rho_mes);
X = H_Mes_Trap_Vec;
N = length(X);

param = inv([N sum(X); sum(X) sum(X.^2)])*[sum(Y);sum(Y.*X)];


H_s = -1/param(2);
Rho_0 = exp(param(1));


Ys = param(2)*X+param(1);


RMS_lin = sqrt(1/N*sum((Ys-Y).^2));
R2_lin = sum((Ys-mean(Y)).^2)/sum((Y-mean(Y)).^2);

Rho_approx = Rho_0*exp(-H_Mes_Trap_Vec/H_s);

V_approx = V_IniM*exp(1/2*B*H_s*(Rho_approx-Rho_approx(1))/sind(-90));

Acc_approx = (1/2.*Rho_approx.*V_approx.^2*S_c*C_do/M_c)';

% RMS de l'approximation au niveau de l'accélération

RMS_V_rel = sqrt(1/N*sum(((V_Mes_Trap_Vec-V_approx)./V_Mes_Trap_Vec).^2));

RMS_acc_abs = sqrt(1/N*sum(((acc_mes-Acc_approx)).^2));

RMS_acc_rel = sqrt(1/N*sum(((acc_mes-Acc_approx)./acc_mes).^2));

%Graph accel
figure()
hold on
title('Accélération en fonction de la hauteur')
plot(H_Mes_Trap_Vec,Acc_approx)
plot(H_Mes_Trap_Vec,acc_mes,'p','color','red')
xlabel('Hauteur en m')
ylabel('Accélération en m/s^2')
legend('Russe','RAA')

%Graph Vitesses
figure()
hold on
title('Vitesse en fonction de la hauteur')
plot(H_Mes_Trap_Vec,V_approx)
plot(H_Mes_Trap_Vec,V_Mes_Trap_Vec,'p','color','red')
xlabel('Hauteur en m')
ylabel('Vitesse en m/s')
legend('Russe','RAA')


%% Loi de guidage
%--------------- A partir de cette ligne en développement-----------------%
%H_FinM/V_FinM1/V_FinM2


%Step 1 - RAA  avec donées de russes

V_RAA = V_IniR*exp(1/2*B*H_s*(Rho_approx-Rho_approx(1))/sin(Gamma_IniR));

% plot(V_RAA,H_Mes_Trap_Vec/1000)
% hold on
% plot(V_Mes_Trap_Vec,H_Mes_Trap_Vec/1000)

RMS_V_RAA = sqrt(1/N*sum(((V_RAA-V_Mes_Trap_Vec)./V_Mes_Trap_Vec).^2));

% 
H = linspace(H_IniM,H_FinM,11000);
Rho_fin = Rho_0*exp(-H_FinM/H_s);
Rho_ini = Rho_0*exp(-H_IniM/H_s);
Rho_RAA = Rho_0*exp(-H/H_s);


%% Section 300

Delta_V_aero300 = V_FinM(2) - sqrt(V_IniM^2+2*Mu_mars*(1/(R_mars+H_FinM)-1/(R_mars+H_IniM)));
Gamma_ref300 = asin(1/2*B*H_s*(Rho_fin-Rho_ini)./(log(1+Delta_V_aero300./V_IniM)));
V_RAA300 = V_IniM*exp(1/2*B*H_s*((Rho_RAA-Rho_RAA(1))./sin(Gamma_ref300)));
Daero_gamma300 = 1/2*Rho_RAA.*V_RAA300.^2*S_c*C_do;

hmin_300 = 10000;
y1_300 = 1/2*Rho_0*S_c*C_do*exp(-hmin_300/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)/sin(Gamma_ref300))-2000;
dy1_300 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmin_300/H_s+B*H_s.*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)./sin(Gamma_ref300))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref300).*exp(-2*hmin_300/H_s+B*H_s*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)./sin(Gamma_ref300)));
it=0;

while abs(y1_300)>1e-8 & it<501
    hmin_300 = hmin_300 - y1_300./dy1_300;
    y1_300 = 1/2*Rho_0*S_c*C_do*exp(-hmin_300/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)/sin(Gamma_ref300))-2000;
    dy1_300 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmin_300/H_s+B*H_s.*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)./sin(Gamma_ref300))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref300).*exp(-2*hmin_300/H_s+B*H_s*(Rho_0*exp(-hmin_300/H_s)-Rho_ini)./sin(Gamma_ref300)));
    it = it+1;
end

hmax_300 = 50000;
y2_300 = 1/2*Rho_0*S_c*C_do*exp(-hmax_300/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)/sin(Gamma_ref300))-2000;
dy2_300 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmax_300/H_s+B*H_s.*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)./sin(Gamma_ref300))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref300).*exp(-2*hmax_300/H_s+B*H_s*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)./sin(Gamma_ref300)));
it=0;

while abs(y2_300)>1e-8 & it<501
    hmax_300 = hmax_300 - y2_300./dy2_300;
    y2_300 = 1/2*Rho_0*S_c*C_do*exp(-hmax_300/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)/sin(Gamma_ref300))-2000;
    dy2_300 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmax_300/H_s+B*H_s.*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)./sin(Gamma_ref300))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref300).*exp(-2*hmax_300/H_s+B*H_s*(Rho_0*exp(-hmax_300/H_s)-Rho_ini)./sin(Gamma_ref300)));
    it = it+1;
end





%% Section 250

Delta_V_aero250 = V_FinM(1) - sqrt(V_IniM^2+2*Mu_mars*(1/(R_mars+H_FinM)-1/(R_mars+H_IniM)));
Gamma_ref250 = asin(1/2*B*H_s*(Rho_fin-Rho_ini)./(log(1+Delta_V_aero250./V_IniM)));
V_RAA250 = V_IniM*exp(1/2*B*H_s*((Rho_RAA-Rho_RAA(1))./sin(Gamma_ref250)));
Daero_gamma250 = 1/2*Rho_RAA.*V_RAA250.^2*S_c*C_do;

hmin_250 = 10000;
y1_250 = 1/2*Rho_0*S_c*C_do*exp(-hmin_250/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)/sin(Gamma_ref250))-2000;
dy1_250 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmin_250/H_s+B*H_s.*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)./sin(Gamma_ref250))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref250).*exp(-2*hmin_250/H_s+B*H_s*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)./sin(Gamma_ref250)));
it=0;

while abs(y1_250)>1e-8 & it<501
    hmin_250 = hmin_250 - y1_250./dy1_250;
    y1_250 = 1/2*Rho_0*S_c*C_do*exp(-hmin_250/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)/sin(Gamma_ref250))-2000;
    dy1_250 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmin_250/H_s+B*H_s.*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)./sin(Gamma_ref250))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref250).*exp(-2*hmin_250/H_s+B*H_s*(Rho_0*exp(-hmin_250/H_s)-Rho_ini)./sin(Gamma_ref250)));
    it = it+1;
end

hmax_250 = 50000;
y2_250 = 1/2*Rho_0*S_c*C_do*exp(-hmax_250/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)/sin(Gamma_ref250))-2000;
dy2_250 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmax_250/H_s+B*H_s.*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)./sin(Gamma_ref250))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref250).*exp(-2*hmax_250/H_s+B*H_s*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)./sin(Gamma_ref250)));
it=0;

while abs(y2_250)>1e-8 & it<501
    hmax_250 = hmax_250 - y2_250./dy2_250;
    y2_250 = 1/2*Rho_0*S_c*C_do*exp(-hmax_250/H_s)*V_IniM^2*exp(B*H_s*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)/sin(Gamma_ref250))-2000;
    dy2_250 = 1/2*S_c*C_do*(-Rho_0*V_IniM^2/H_s.*exp(-hmax_250/H_s+B*H_s.*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)./sin(Gamma_ref250))-V_IniM^2*B*Rho_0^2./sin(Gamma_ref250).*exp(-2*hmax_250/H_s+B*H_s*(Rho_0*exp(-hmax_250/H_s)-Rho_ini)./sin(Gamma_ref250)));
    it = it+1;
end



%%


Rho_min_300 = Rho_0*exp(-hmin_300/H_s);
Rho_max_300 = Rho_0*exp(-hmax_300/H_s);
Rho_min_250 = Rho_0*exp(-hmin_250/H_s);
Rho_max_250 = Rho_0*exp(-hmax_250/H_s);


V_min300 = sqrt(2*2000/S_c/C_do/Rho_min_300);
V_max300 = sqrt(2*2000/S_c/C_do/Rho_max_300);
V_min250 = sqrt(2*2000/S_c/C_do/Rho_min_250);
V_max250 = sqrt(2*2000/S_c/C_do/Rho_max_250);


V_moy300 = 1/2*(V_min300+V_max300);
V_moy250 = 1/2*(V_min250+V_max250);


DT_300 = (hmin_300-hmax_300)/sin(Gamma_ref300)/V_moy300
DT_250 = (hmin_250-hmax_250)/sin(Gamma_ref250)/V_moy250

figure()
hold on
plot(H,Daero_gamma300)
xline(hmin_300,':','color','red');
xline(hmax_300,':','color','red');
yline(2000,':','color','red');

figure()
hold on
plot(H,Daero_gamma250)
xline(hmin_250,':','color','red');
xline(hmax_250,':','color','red');
yline(2000,':','color','red');


disp('end')

save App6Problematique_Identification.mat Rho_0 H_s
