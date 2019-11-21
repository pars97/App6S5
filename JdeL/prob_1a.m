% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (révision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (révision)
% 
% Problème no 1: Intégrateur numérique,  stabilité et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un intégrateur numérique sur MATLAB
% (b) connaître les méthodes d'intégration explicite et prédiction-correction
% (c) savoir reconnaître une instabilité numérique et en trouver la cause
% (d) connaître et savoir utiliser le lien entre l'erreur de discrétisation 
%     et le pas d'intégration.

% Problème no 1a : Équation différentielle marginalement stable qui est déstabilisée
% 				   par les erreurs de discrétisation.
%
% Dans la partie 1a, on trouve la solution numérique d'un problème dynamique
% avec 3 différents intégrateurs numériques:
% (1) Euler explicite
% (2) Euler implicite
% (3) Runge-Kutta
% On compare ensuite leur stabilité et précision en fonction du pas d'intégration.

% NOTE: Les intégrateurs numériques de type Euler sont rarement utilisés en
% pratique. Ils ne sont utilisés ici que pour des raisons pédagogiques: leur
% mauvaise performance permet d'illuster plus facilement des problèmes de
% stabilité et de convergence parfois rencontrés avec des intégrateurs
% numériques plus performants quand les équations différentielles à
% intégrer comportent des problèmes particuler ('stiffness' par exemple). 
% --------------------------------------------------------------------------------
% 
  clc
  close all  
  clear 

  figure(1)		% Commence nouvelle figure
  N = 10;		% Nombre de différents pas d'intégration
% dt = 2.0;		% Pas d'intégration 
  dt = 1.0;		% Pas d'intégration 
  y0 = [1 0];	% Conditions initiales
% y0 = [1];	    % Conditions initiales
  
% Boucle pour calculer les résultats numériques pour chaque pas d'intégration
  
  for n=1:N,
     
%    Intégration avec méthode de Euler     
     
     pas(n) = dt;									% pas d'intégration à utiliser
     [t, yex]=eulerex('eqn1a',0, dt, 4.5, y0);		% intégration avec Euler explicite
     disp([n, max(t)])
     [t, ypc]=eulerpc('eqn1a',0, dt, 4.5, y0);		% intégration avec Euler prédiction-correction
     disp([n, max(t)])

%    Integration avec méthode Runge-Kutta 
%    (On utilise une grande tolérance pour avoir un pas fixe égal à dt)
     
     options = odeset('abstol',1000,'reltol',1000,'maxstep',dt);
     [t, yrk] = ode45('eqn1a',t,y0,options);
     
%    Solution exacte et indice de la valeur finale f
     
     y = cos(3*t);
%    y = y0*exp(2*t);
     f = length(t);
     
%    Graphique des erreurs de position et calcul de l'erreur 

     subplot(3,1,1)
     title('Comparaison: solutions numérique et exacte')
     plot(t, [yex(:,1), y])
     ylabel('Y EX')
     hold on
     errex(n) = abs(yex(f,1)-y(f));
%    errex(n) = sqrt((yex(:,1)-y)'*(yex(:,1)-y));
     
     subplot(3,1,2)
     plot(t, [ypc(:,1), y])
     ylabel('Y PC')
     hold on
     errpc(n) = abs(ypc(f,1)-y(f));
%    errpc(n) = sqrt((ypc(:,1)-y)'*(ypc(:,1)-y));
     
     subplot(3,1,3)
     plot(t, [yrk(:,1), y])
     ylabel('Y RK')
     hold on
     errrk(n) = abs(yrk(f,1)-y(f));
%    errrk(n) = sqrt((yrk(:,1)-y)'*(yrk(:,1)-y));
     
%    Réduction du pas d'intégration pour la prochaine boucle
       
     dt = dt/2;

  end

  figure(1)
  subplot(3,1,1), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  subplot(3,1,2), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  subplot(3,1,3), xlabel('Temps (s)'), grid, axis([0.0 4.5 -2 2])
  hold off
  
% Graphique des erreurs sur 2eme figure
  
  figure(2)

  subplot(3,1,1)
  plot(pas',errex') 
  title('Erreur en fonction du pas d''intégration')
  ylabel('EULER EX')
  grid on
  
  subplot(3,1,2)
  plot(pas',errpc')
  ylabel('EULER PC')
  grid on
  
  subplot(3,1,3)
  plot(pas',errrk')
  ylabel('RK')
  xlabel('Pas d''intégration')
  grid on

% Calcul de l'ordre des méthodes avec les "vraies" erreurs, utilisant la solution analytique exacte  
  
% ordre_ex = log(errex(2:N)./errex(1:N-1))./log(pas(2:N)./pas(1:N-1));
% ordre_pc = log(errpc(2:N)./errpc(1:N-1))./log(pas(2:N)./pas(1:N-1));
% ordre_rk = log(errrk(2:N)./errrk(1:N-1))./log(pas(2:N)./pas(1:N-1));

  ordre_ex = log(errex(1:N-1)./errex(N))./log(pas(1:N-1)./pas(N));
  ordre_pc = log(errpc(1:N-1)./errpc(N))./log(pas(1:N-1)./pas(N));
  ordre_rk = log(errrk(1:N-1)./errrk(N))./log(pas(1:N-1)./pas(N));

%   ordre_ex = log(errex(1:N-1)./errex(2:N))./log(pas(1:N-1)./pas(2:N));
%   ordre_pc = log(errpc(1:N-1)./errpc(2:N))./log(pas(1:N-1)./pas(2:N));
%   ordre_rk = log(errrk(1:N-1)./errrk(2:N))./log(pas(1:N-1)./pas(2:N));
  
  figure(3)
  subplot(3,1,1)
  plot(pas(1:N-1)',ordre_ex')
  ylabel('EULER EX')
  title('Calcul de l''ordre de l''intégrateur: détection de convergence')
  axis([0.0 0.8 0 2])
  grid on
  
  subplot(3,1,2)
  plot(pas(1:N-1)',ordre_pc')
  ylabel('EULER PC')
  axis([0.0 0.8 0 3])
  grid on
  
  subplot(3,1,3)
  plot(pas(1:N-1)',ordre_rk')
  ylabel('ODE45')
  xlabel('Pas d''intégration')
  axis([0.0 0.8 4.5 5.5])
  grid on

  figure(4)
  subplot(3,1,1)
  plot(log(pas(1:N-1)'/pas(N)),log(errex(1:N-1)'/errex(N)))
  ylabel('EULER EX')
  title('Erreur (rapport logarithmique)')
  grid on
  
  subplot(3,1,2)
  plot(log(pas(1:N-1)'/pas(N)),log(errpc(1:N-1)'/errex(N)))
  ylabel('EULER PC')
  grid on
  
  subplot(3,1,3)
  plot(log(pas(1:N-1)'/pas(N)),log(errrk(1:N-1)'/errex(N)))
  ylabel('ODE45')
  xlabel('Pas d''intégration (rapport logarithmique)')
  grid on


