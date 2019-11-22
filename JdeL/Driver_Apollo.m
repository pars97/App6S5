% S4-APP4 JdeLafontaine 17 juin 2008
% S5-APP6 JdeLafontaine 16 novembre 2019 (révision)
% 
% Intégrateur numérique Runge Kutta

% (1) Effet de la 'reltol' de l'intégrateur sur le nombre de pas et la
%     précision de la trajectoire (qui doit se refermer sur elle-même)
% (2) Variation du pas d'intégration selon l'amplitude des derivées (près de la Terre et de la Lune)
% (3) Utilisation de ODE45

% Etat: z = [x, y, dx/dt, dy/dt]

  clc
  clear
  close all

% Conditions initiales et temps final

  u = 1.0/(82.45);
  z0 = [1.2, 0.0, 0.0, -1.04935751];
  tspan = [0, 6.19];
  
% Place la Lune et la Terre sur le graphique

  plot((1-u),0,'rp', 'Markersize',5)
  hold on, plot(-u,0,'ro', 'Markersize',5)

% 1ère intégration avec paramètres nominaux

  reltol1 = 1e-03;
  options = odeset('abstol' ,1e-06, 'reltol', reltol1);
  [t, z] = ode45('apollo', tspan, z0, options);

  plot(z(:,1), z(:,2),'.b', 'Markersize',10)
  Lz = length(z);
  plot(z(Lz,1), z(Lz, 2),'xb', 'Markersize',10)
  
% 2ème intégration avec paramètres reltol = 1e-06

  reltol2 = 1e-06;
  options = odeset('abstol', 1e-06, 'reltol', reltol2);
  [t, z] = ode45('apollo', tspan, z0, options);
  
  plot(z(:,1), z(:,2),'.g', 'Markersize',10)
  Lz = length(z);
  plot(z(Lz,1), z(Lz, 2),'xg', 'Markersize',10)
  
  xlabel('Position X', 'Fontsize',15)
  ylabel('Position Y', 'Fontsize',15)
  title('Capsule Apollo', 'Fontsize',15)

  legend('Lune', 'Terre', ['Précision: ', num2str(reltol1)], 'Fin trajectoire', ['Précision: ', num2str(reltol2)], 'Fin trajectoire')
