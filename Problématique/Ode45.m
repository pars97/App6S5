
% Etat: z = [v,gamma,h,s,theta,q]

  clc
  clear
  close all
  

% Conditions initiales et temps final
  z0 = [6100 deg2rad(-20.5) 120000 0 deg2rad(-80) 0];
  tspan = [0,150];

% 1ère intégration avec paramètres nominaux

  reltol1 = 1e-02;
  options = odeset('abstol' ,1e-02, 'reltol', reltol1);
  [t, z] = ode45('ass', tspan, z0, options);
hold on
plot(t,z(:,1))
plot(t,z(:,3))
yline(10000);
legend('Vitesse','Hauteur')
