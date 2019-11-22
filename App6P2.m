clear all 
close all
clc
%% APP6P2

e = [0.0843 0.02603 0.01048 0.00319 0.00040];
delta_t = [0.05 0.04 0.03 0.02 0.01];


p1 = log(e(1:end-1)./e(end))./log(delta_t(1:end-1)./delta_t(end));
p2 = log(e(1:end-1)./e(2:end))./log(delta_t(1:end-1)./delta_t(2:end));


plot(delta_t(1:end-1),p1)
plot(delta_t(1:end-1),p2)






