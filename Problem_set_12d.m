close all
clear all

%Declaring variables
E0 = 1e4; %[MPa]
E1 = 1e4; %[MPa]
eta1 = 1e7; %[MPa/s]
tau1 = eta1/E1; %[1/s] - relaxation time
t0 = 0.1*tau1; %(relaxation time)/10 gives 10 measurements before relaxation
dt = 0.1*t0;
t = 0:dt:21*t0;
sigma0 = E0/100;
epsilon0 = sigma0/E0; %epsilon0 = 1/100 as epsilon0 > 0
N = length(t); %number of measurement points

%creating empty vectors to save history
sigma_1 = zeros(1,N);
sigma_2 = zeros(1,N);

epsilon_1 = zeros(1,N);
epsilon_1i = zeros(1,N);
epsilon_2 = zeros(1,N);
epsilon_2i = zeros(1,N);

