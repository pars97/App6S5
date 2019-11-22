%% NO14
% Etat: z = [x, y, dx/dt, dy/dt]

  clc
  clear
  close all

% Conditions initiales et temps final

  u = 0;
  z0 = [-3 7.832];
  tspan = [0, 10];


  reltol1 = 1e-08;
  options = odeset('abstol' ,1e-08, 'reltol', reltol1);
  [t, z] = ode45('no14', tspan, z0, options);

  plot(z(:,1),z(:,2))
 








