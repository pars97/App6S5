% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (révision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (révision)
% 
% Probleme no 1: Intégrateur numérique,  stabilité et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un intégrateur numérique sur MATLAB
% (b) connaître les méthodes d'intégration explicite et prédiction-correction
% (c) savoir reconnaître une instabilité numérique et en trouver la cause
% (d) connaître et savoir utiliser le lien entre l'erreur de discrétisation 
%     et le pas d'integration.

% Problème no 1c : Équation différentielle instable.
%
% Dans cette partie 1c, on trouve la solution numérique d'un problème dynamique de 
% nature instable. On utilise 3 différents pas d'intégration. 

% NOTE: Les intégrateurs numériques de type Euler sont rarement utilisés en
% pratique. Ils ne sont utilisés ici que pour des raisons pédagigiques: leur
% mauvaise performance permet d'illuster plus facilement des problèmes de
% stabilité et de convergence parfois rencontrés avec des intégrateurs
% numériques plus performants quand les équations différentielles à
% intégrer comportent des problèmes particuler ('stiness' par exemple). 

  clc
  close all
  clear
  
  pas = [0.002, 0.001, 0.0005];
  for n=1:3,
     [t, yex]=eulerex('eqn1c',0, pas(n), 3.0 ,0.08);
     [t, ypc]=eulerpc('eqn1c',0, pas(n), 3.0 ,0.08);
     
     y = t.^2 + 0.4*t + 0.08;										% Solution exacte

     figure(1)							
     subplot(3,2,2*n-1)
     plot(t, [yex, y])
     if(n==1), title('EULEREX'), end
     ylabel(['Pas = ', num2str(pas(n))])

     subplot(3,2,2*n)
     plot(t, [ypc, y])
     if (n==1), title('EULERPC'), end
     ylabel(['Pas = ', num2str(pas(n))])
  end
  figure(1), subplot(3,2,1), legend('EulerEX', 'Exacte')
  figure(1), subplot(3,2,2), legend('EulerPC', 'Exacte')
  figure(1), subplot(3,2,5), xlabel('Temps (s)'), subplot(3,2,6), xlabel('Temps (s)')




