close all
clear all
clc

% Graph 250
load ('Ode45Donnee')

plot(t,[z(:,1)-z1(:,1), z(:,3)]-z1(:,3))
xline(131.2)
legend('Vitesse','Hauteur')
title ('V_fin = 250m/s')







clear all

% Graph 300
load ('Ode45Donneesvf300')
plot(t,[z(:,1), z(:,3)])
xline(114.2)
legend('Vitesse','Hauteur')
title ('V_fin = 300m/s')
clear all





