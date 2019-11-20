clear all 
close all
clc
%% Approximation à 2 paramètres
 
%E5c)

x = [1 2 3 4 5];
y = [11.4 6.1 3.6 3.3 2.1];

Y = y;
X = 1./x;

N = length(X);

param = inv([N,sum(X); sum(X),sum(X.^2)])*[sum(Y);sum(X.*Y)];

alpha = param(1);
beta = param(2);


gx = alpha+beta./x;

erreur_quad = sum((gx-y).^2);




