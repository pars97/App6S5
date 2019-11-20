close all 
clear all
clc

%% numero 8
% Fonction plot 
x = linspace (0,5,1000);
y = x.^3-6*x.^2 +7*x+2;
dy = 3*x.^2-12*x+7;


%% Valeur de départ de 0.85
it = 0;
Xn = 0.85;
y = Xn^3-6*Xn^2 +7*Xn+2;
dy = 3*Xn^2-12*Xn+7;


while abs(y)>1e-8 && it<501
    Xn = Xn - y/dy;
    y = Xn^3-6*Xn^2 +7*Xn+2;
    dy = 3*Xn^2-12*Xn+7;
    it = it+1;
end

disp(Xn)
disp(it)




%% Valeur de départ de 1.10

it = 0;
Xn = 1.10;
y = Xn^3-6*Xn^2 +7*Xn+2;
dy = 3*Xn^2-12*Xn+7;
while abs(y)>1e-8 && it<501
    Xn = Xn - y/dy;
    y = Xn^3-6*Xn^2 +7*Xn+2;
    dy = 3*Xn^2-12*Xn+7;
    it = it+1;
end

disp(Xn)
disp(it)



%% Valeur de départ 0.709

it = 0;
Xn = 0.85;
y = Xn^3-6*Xn^2 +7*Xn+2;
dy = 3*Xn^2-12*Xn+7;
while abs(y)>1e-8 && it<501
    Xn = Xn - y/dy;
    y = Xn^3-6*Xn^2 +7*Xn+2;
    dy = 3*Xn^2-12*Xn+7;
    it = it+1;
end

disp(Xn)
disp(it)



%% Valeur de départ de 1
it = 0;
Xn = 1;
y = Xn^3-6*Xn^2 +7*Xn+2;
dy = 3*Xn^2-12*Xn+7;
while abs(y)>1e-8 && it<501
    Xn = Xn - y/dy;
    y = Xn^3-6*Xn^2 +7*Xn+2;
    dy = 3*Xn^2-12*Xn+7;
    it = it+1;
end

disp(Xn)
disp(it)






