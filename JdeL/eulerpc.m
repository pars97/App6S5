  function [t,y]=eulerpc(fct,t0,dt,tf,y0)
%
% Fonction qui intègre les équations différentielles "fct" selon la
% méthode d'Euler prédiction-correction, de t0 a tf avec pas dt.
%
% Utilisation: [t,y]=eulerpc(fct,y0,t0,dt,tf)
%
% Entrées:
% fct : variable alpha-numérique avec le nom de la fonction f(y,t) à intégrer
% y0  : condition initiale
% t0  : temps initial
% dt  : pas d'intégration (fixe)
% tf  : temps final
%
% Sorties:
% y   : solution numérique avec une rangée par pas d'intégration et 
%       une colonne par variables de sortie (nbre de col = ordre du système).
% t   : matrice-colonne temps de t0 à tf par intervals dt
%--------------------------------------------------------------------
%
  t = [t0:dt:tf]';		% Calcul de la mat-col temps
  nstep=size(t)-1 ;     % Calcul du nombre de pas
  y(1,:)=y0(:)';			% Initialisation
  for n=1:nstep,			% Boucle pour propager l'état du système
    fn = feval(fct,t(n),y(n,:))';
    yp = y(n,:)+dt*fn;			% Prédiction
    fp = feval(fct,t(n+1),yp)';
    y(n+1,:)=y(n,:)+dt*(fn+fp)/2;	% Correction
  end
%








