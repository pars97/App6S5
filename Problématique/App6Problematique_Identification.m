close all
clear all
clc
addpath ../
load('Accelero_Data_from_Moscow');
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
E_H_Mes_Trap = H_V_Trap^2/12*(DpVmes_b-DpVmes_a)

%intégrale point par point
H_Mes_Trap_Vec(1)=H_IniR;


for n=2:1:length(Temps)
    H_Mes_Trap_Vec(n) = H_IniR - H_V_Trap/2*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(n)+2*sum(V_Mes_Trap_Vec(2:n-1)));
end

% Détermination de la vitesse selon l'accélération par méthode de Simpson
DtHmes_a = (V_Mes_Trap_Vec(4)-3*V_Mes_Trap_Vec(3)+3*V_Mes_Trap_Vec(2)-V_Mes_Trap_Vec(1))/H_V_Trap^3;
DtHmes_b = (V_Mes_Trap_Vec(end)-3*V_Mes_Trap_Vec(end-1)+3*V_Mes_Trap_Vec(end-2)-V_Mes_Trap_Vec(end-3))/H_V_Trap^3;

H_Mes_Simp = H_IniR + H_V_Trap/3*(V_Mes_Trap_Vec(1)+V_Mes_Trap_Vec(end)+4*sum(V_Mes_Trap_Vec(2:2:end-1))+2*sum(V_Mes_Trap_Vec(3:2:end-1)));
E_H_Mes_Simp = -H_V_Trap^4/180*(V_Mes_Trap_Vec(end)-V_Mes_Trap_Vec(1))*(DtHmes_b-DtHmes_a)

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

Acc_approx = (1/2.*Rho_approx.*V_Mes_Trap_Vec.^2*S_c*C_do/M_c)';

% RMS de l'approximation au niveau de l'accélération

RMS_acc_rel = sqrt(1/N*sum(((acc_mes-Acc_approx)).^2))

RMS_acc_rel = sqrt(1/N*sum(((acc_mes-Acc_approx)./acc_mes).^2))

% plot(H_Mes_Trap_Vec,[acc_mes,Acc_approx])

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



Delta_V_aero = V_FinM - sqrt(V_IniM^2+2*Mu_mars*(1/(R_mars+H_FinM)-1/(R_mars+H_IniM)));

Gamma_ref = asin(1/2*B*H_s*(Rho_fin-Rho_ini)./(log(1+Delta_V_aero./V_IniM)));

Rho_RAA = Rho_0*exp(-H/H_s);

V_RAA = V_IniM*exp(1/2*B*H_s*((Rho_RAA-Rho_RAA(1))./sin(Gamma_ref)));

Daero_gamma = 1/2*Rho_RAA.*V_RAA.^2*S_c*C_do;


%% Newton raphson
it1 = 0;
h1 = H(10000);
it2 = 0;
h2 = H(5000);

y1 = 1/2*Rho_0*exp(-h1/H_s)*S_c*C_do*(V_IniM*exp(1/2*B*H_s*((Rho_0*exp(-h1/H_s)-Rho_0*exp(-H(1)/H_s))./sin(Gamma_ref)))).^2-2000;
dy1 = - (6*C_do*Rho_0^2*S_c*V_IniM^2*exp(-(12*H_s*(Rho_0*exp(-H(1)/H_s) - Rho_0*exp(-h1/H_s)))./(625*sin(Gamma_ref)))*exp(-(2*h1)/H_s))./(625*sin(Gamma_ref)) - (C_do*Rho_0*S_c*V_IniM^2*exp(-(12*H_s*(Rho_0*exp(-H(1)/H_s) - Rho_0*exp(-h1/H_s)))./(625*sin(Gamma_ref)))*exp(-h1/H_s))/(2*H_s);

while abs(y1)>1e-8 & it1<500
    h1 = h1 - y1./dy1;
    y1 = 1/2*Rho_0*exp(-h1/H_s)*S_c*C_do.*(V_IniM*exp(1/2*B*H_s*((Rho_0*exp(-h1/H_s)-Rho_0*exp(-H(1)/H_s))./sin(Gamma_ref)))).^2-2000;
    dy1 = - (6*C_do*Rho_0^2*S_c*V_IniM^2.*exp(-(12*H_s*(Rho_0.*exp(-H(1)/H_s) - Rho_0.*exp(-h1/H_s)))./(625*sin(Gamma_ref))).*exp(-(2*h1)/H_s))./(625*sin(Gamma_ref)) - (C_do*Rho_0*S_c*V_IniM^2.*exp(-(12*H_s*(Rho_0.*exp(-H(1)/H_s) - Rho_0.*exp(-h1/H_s)))./(625*sin(Gamma_ref))).*exp(-h1/H_s))/(2*H_s);
    it1 = it1+1;
end

y2 = 1/2*Rho_0*exp(-h2/H_s)*S_c*C_do*(V_IniM*exp(1/2*B*H_s*((Rho_0*exp(-h2/H_s)-Rho_0*exp(-H(1)/H_s))./sin(Gamma_ref)))).^2-2000;
dy2 = - (6*C_do*Rho_0^2*S_c*V_IniM^2*exp(-(12*H_s*(Rho_0*exp(-H(1)/H_s) - Rho_0*exp(-h2/H_s)))./(625*sin(Gamma_ref)))*exp(-(2*h2)/H_s))./(625*sin(Gamma_ref)) - (C_do*Rho_0*S_c*V_IniM^2*exp(-(12*H_s*(Rho_0*exp(-H(1)/H_s) - Rho_0*exp(-h2/H_s)))./(625*sin(Gamma_ref)))*exp(-h2/H_s))/(2*H_s);

while abs(y2)>1e-8 & it2<500
    h2 = h2 - y2./dy2;
    y2 = 1/2*Rho_0*exp(-h2/H_s)*S_c*C_do.*(V_IniM*exp(1/2*B*H_s*((Rho_0*exp(-h2/H_s)-Rho_0*exp(-H(1)/H_s))./sin(Gamma_ref)))).^2-2000;
    dy2 = - (6*C_do*Rho_0^2*S_c*V_IniM^2.*exp(-(12*H_s*(Rho_0.*exp(-H(1)/H_s) - Rho_0.*exp(-h2/H_s)))./(625*sin(Gamma_ref))).*exp(-(2*h2)/H_s))./(625*sin(Gamma_ref)) - (C_do*Rho_0*S_c*V_IniM^2.*exp(-(12*H_s*(Rho_0.*exp(-H(1)/H_s) - Rho_0.*exp(-h2/H_s)))./(625*sin(Gamma_ref))).*exp(-h2/H_s))/(2*H_s);
    it2 = it2+1;
end

figure()
hold on
plot(H,Daero_gamma(1,:))
xline(h1(1),':','color','red');
xline(h2(1),':','color','red');
yline(2000,':','color','red');

figure()
hold on
plot(H,Daero_gamma(2,:))
xline(h1(2),':','color','red');
xline(h2(2),':','color','red');
yline(2000,':','color','red');

%%

Hdot = sqrt(2*Daero_gamma./(Rho_0.*exp(-H/H_s).*S_c*C_do)).*sin(Gamma_ref);
Hdot1_moy = mean(Hdot(1,(h1(1)/10:h2(1)/10)));
Hdot2_moy = mean(Hdot(1,(h1(2)/10:h2(2)/10)));
Dt_1 = (h1(1)-h2(1))/Hdot1_moy;
Dt_2 = (h1(2)-h2(2))/Hdot2_moy;

disp('end')


