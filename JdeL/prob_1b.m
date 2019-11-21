% S4-APP4 JdeLafontaine 17 juin 2008
% S4-APP6 JdeLafontaine 18 juillet 2013 (révision)
% S5-APP6 JdeLafontaine 16 novembre 2019 (révision)

% Problème no 1: Intégrateur numérique,  stabilité et convergence
% 
% Objectifs: 
% (a) savoir programmer et utiliser un intégrateur numérique sur MATLAB
% (b) connaître les méthodes d'intégration explicite et prédiction-correction
% (c) savoir reconnaître une instabilité numérique et en trouver la cause
% (d) connaître et savoir utiliser le lien entre l'erreur de discrétisation 
%     et le pas d'integration.

% Problème no 1b : Équation différentielle stable mais "stiff" (ses deux valeurs propres
%                  sont très différentes). L'intégrateur numérique est 
%                  conditionnellement stable (dépend du pas d'intégration).
%
% Dans cette partie 1b, on trouve la solution numérique d'un problème dynamique de 
% nature "stiff" i.e. avec des valeurs propres très différentes. On utilise 5 différents
% pas d'intégration et on observe ensuite la stabilité en fonction du pas d'intégration.

% À noter que pour le pas = 0.002, la stabilité est marginale (pôles sur cercle unitaire)

% NOTE: Les intégrateurs numériques de type Euler sont rarement utilisés en
% pratique. Ils ne sont utilisés ici que pour des raisons pédagigiques: leur
% mauvaise performance permet d'illuster plus facilement des problèmes de
% stabilité et de convergence parfois rencontrés avec des intégrateurs
% numériques plus performants quand les équations différentielles à
% intégrer comportent des problèmes particuler ('stiness' par exemple). 
  
  clear errex errpc
  close all
  clc

  pas = [0.010, 0.005, 0.002, 0.0015, 0.0010];
  for n=1:5,
     [t, yex]=eulerex('eqn1b',0, pas(n), 0.05,[1 1]);
     [t, ypc]=eulerpc('eqn1b',0, pas(n), 0.05,[1 1]);
     y1 = 4*exp(-t) - 3*exp(-1000*t);
     y2 =-2*exp(-t) + 3*exp(-1000*t);

     figure(1)							% Sortie 1
     subplot(5,2,2*n-1)
     plot(t, [yex(:,1), y1])
     ylabel(['Pas = ', num2str(pas(n))])
     if(n==1), title('EULEREX'), end
     subplot(5,2,2*n)
     plot(t, [ypc(:,1), y1])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULERPC'), end
   
     figure(2)							% Sortie 2
     subplot(5,2,2*n-1)
     plot(t, [yex(:,2), y2])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULEREX'), end
     subplot(5,2,2*n)
     plot(t, [ypc(:,2), y2])
     ylabel(['Pas = ', num2str(pas(n))])
     if (n==1), title('EULERPC'), end
 
  end

  figure(1), subplot(5,2,1), legend('EulerEX', 'exacte')
  figure(1), subplot(5,2,2), legend('EulerPC', 'exacte')

  figure(2), subplot(5,2,1), legend('EulerEX', 'exacte')
  figure(2), subplot(5,2,2), legend('EulerPC', 'exacte')



