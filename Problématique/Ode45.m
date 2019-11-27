
% Etat: z = [v,gamma,h,s,theta,q]

  clc
  clear
  close all
  
load('App6Problematique_Identification')
run('constantes')
% Conditions initiales et temps final
  z0 = [6100 deg2rad(-20.5) 120000 0 deg2rad(-80) 0 0];
  tspan = [0,130];

% 1ère intégration avec paramètres nominaux

  reltol1 = 1e-08;
  options = odeset('abstol' ,1e-08, 'reltol', reltol1);
  [t, z] = ode45('ass', tspan, z0, options);
  %[t1, z1] = ode45('ass1', tspan, z0, options);
  %%

  
P_dyn  = (1/2)*Rho_0*exp(-z(:,3)/H_s).*z(:,1).^2;
D_aero = P_dyn.*S_c*C_do;
L_aero = P_dyn.*S_c*C_Lalpha.*(z(:,5)-z(:,2));

r_fin = R_mars + H_FinM;
DVA = V_FinM2 - sqrt(z(:,1).^2 + ((2*Mu_mars).*((1/r_fin - 1./(R_mars + z(:,3))))));

Rho_fin = Rho_0*exp(-H_FinM/H_s);
Rho = Rho_0*exp(-z(:,3)/H_s);
Gamma_ref = asin((1/2)*B*H_s*((Rho_fin - Rho) ./ log(1 + DVA./z(:,1))));
%%
  close all
figure()
title('Angle gamma en fonction du temps (300m/s)')
hold on
plot(t,rad2deg(z(:,2)))
% plot(t,rad2deg(Gamma_ref))
xlabel('Temps(s)')
ylabel('Angle gamma(°)')
% axis([0 t(end) -25 -10])
legend('Gamma','Gamma Ref')

figure()
title('Angle gamma en fonction de la hauteur (300m/s)')
hold on
plot(z(:,3)/1000,rad2deg(z(:,2)))
xlabel('Hauteur(km)')
ylabel('Angle gamma(°)')

figure()
title('Vitesse en fonction du temps (300m/s)')
hold on
plot(t,(z(:,1)))
xlabel('Temps(s)')
ylabel('Vitesse(m/s)')


figure()
title('Vitesse en fonction de la hauteur (300m/s)')
hold on
plot(z(:,3)/1000,z(:,1))
xlabel('Hauteur(km)')
ylabel('Vitesse(m/s)')
xline(10);

figure()
title('Hauteur en fonction du temps (300m/s)')
hold on
plot(t,z(:,3)/1000)
ylabel('Hauteur(km)')
xlabel('Temps(s)')

figure()
title('Angle theta et alpha en fonction du temps (300m/s)')
hold on
plot(t,rad2deg(z(:,5))) 
plot(t,rad2deg(z(:,5)-z(:,2)))
xlabel('Temps(s)')
ylabel('Angle(°)')
legend('Angle theta','Angle alpha')

figure()
title('Vitesse angulaire theta en fonction du temps (300m/s)')
hold on
plot(t,rad2deg(z(:,6)))
% axis([0 t(end) -5 5])
xlabel('Temps(s)')
ylabel('Vitesse angulaire(°/s)')


position = find(z(:,7)>0,1, 'first');

figure()
title('Daero et Laero en fonction du temps (300m/s)')
hold on
plot(t,D_aero)
plot(t,L_aero)
xlabel('Temps(s)')
ylabel('Force sur la capsule (N)')
xline(t(position),':','color','red');
xline(t(position)+max(z(:,7)),':','color','red');
yline(2000,'--','color','red');
legend('Daero','Laero')

save Ode300A 


