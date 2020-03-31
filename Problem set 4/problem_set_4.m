%Problem set 4: Tensile test
clear all 
close all
%task 1- checking the scatter in the measured force-displacement curves:
%importing test data from test 1, force F [N] vs. displacement [mm]
T1 = readtable('test-data.xlsx','Range','A19:D1508');
T2 = readtable('test-data.xlsx','Range','F19:I1508');
T3 = readtable('test-data.xlsx','Range','L19:O1508');

%plotting the results of cross-head displacement
figure();
hold on;
title('Cross-head displacement');
xlabel('Displacement [mm]');
ylabel('Force [N]');
xlim([0 30]);
grid on;
plot(T1{:,3},T1{:,2});
plot(T2{:,3},T2{:,2});
plot(T3{:,3},T3{:,2});
leg1 = legend("Test 1", "Test 2", "Test 3");
set(leg1, 'Location', 'SouthEast')
saveas(gcf,'task1_cross_head_disp.png');
hold off;

%plotting the results for extensometer
figure();
hold on;
title('Extensometer');
xlabel('Displacement [mm]');
ylabel('Force [N]');
xlim([0 5]);
plot(T1{:,4},T1{:,2});
plot(T2{:,4},T2{:,2});
plot(T3{:,4},T3{:,2});
leg2 = legend("Test 1", "Test 2", "Test 3");
set(leg2, 'Location', 'best')
saveas(gcf,'task1_extensometer.png');
hold off;

%% Task 2- stress-strain relationship

%declares variables
L0_s = 70; %initial length of specimen [mm]
L0_e = 40; %initial length of extensometer [mm]
t = 3; %specimen initial thichness [mm]
h = 12.5; %specimen initial height [mm]
A0 = h*t; %specimen initial cross sectional area [mm^2], assumes same area for extensometer and specimen

%creates the strain vector for the tests of the specimen:
%Formula: epsilon = (L0-L)/L0 = displacement/initial length
eps1_s = T1{:,3}/L0_s; %strain of specimen, test 1 
eps2_s = T2{:,3}/L0_s; %strain of specimen, test 2 
eps3_s = T3{:,3}/L0_s; %strain of specimen, test 3

eps_s = [eps1_s;eps2_s;eps3_s];

%creates the strain vector for the tests of the extensometer:
%Formula: epsilon = (L0-L)/L0 = displacement/initial length
eps1_e = T1{:,4}/L0_e; %strain of extensometer, test 1 
eps2_e = T2{:,4}/L0_e; %strain of extensometer, test 2 
eps3_e = T3{:,4}/L0_e; %strain of extensometer, test 3

%creates the stress vector for the tests:
%Formula: sigma = F/A0 = force/initial area
sigma1 = T1{:,2}/A0;
sigma2 = T2{:,2}/A0;
sigma3 = T3{:,2}/A0;

sigma = [sigma1;sigma2;sigma3];

%Plots the stress-strain relationship of test 1:
figure();
hold on;
grid on;
plot(eps1_s,sigma1);
plot(eps1_e,sigma1);
title('Stress-strain relationship, test 1');
xlabel('Engineering strain');
ylabel('Engineering stress, [MPa]');
leg3 = legend("Cross-head", "Extensometer");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task2_stress_strain_1.png');
hold off;

%Plots the stress-strain relationship of test 2:
figure();
hold on;
grid on;
plot(eps2_s,sigma2);
plot(eps2_e,sigma2);
title('Stress-strain relationship, test 2');
xlabel('Engineering strain');
ylabel('Engineering stress, [MPa]');
leg3 = legend("Cross-head", "Extensometer");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task2_stress_strain_2.png');
hold off;

%Plots the stress-strain relationship of test 3:
figure();
hold on;
grid on;
plot(eps3_s,sigma3);
plot(eps3_e,sigma3);
title('Stress-strain relationship, test 3');
xlabel('Engineering strain');
ylabel('Engineering stress, [MPa]');
leg3 = legend("Cross-head", "Extensometer");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task2_stress_strain_3.png');
hold off;

%Plots all tests in the same plot:
figure();
hold on;
grid on;
plot(eps1_s,sigma1);
plot(eps1_e,sigma1);
plot(eps3_s,sigma3);
plot(eps3_e,sigma3);
plot(eps2_s,sigma2);
plot(eps2_e,sigma2);
title('Stress-strain relationship, all tests');
xlabel('Engineering strain');
ylabel('Engineering stress, [MPa]');
leg3 = legend("Cross-head_1", "Extensometer_1","Cross-head_2", "Extensometer_2", "Cross-head_3", "Extensometer_3");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task2_stress_strain_all.png');
hold off;

%% Task 3-Correction of results
%variables
L0_s2 = 80; %[mm]

%creates the updated strain vector for the tests of the specimen:
%Formula: epsilon = (L0-L)/L0 = displacement/initial length
eps1_s2 = T1{:,3}/L0_s2; %strain of specimen, test 1 

%The strain from the extensometer reimains the same:
%eps1_e

%creates the stress vector for the tests:
%Formula: sigma = F/A0 = force/initial area
sigma_ce = T1{:,2}/A0;

%finds delta epsilon- deviation of start point of extensometer and specimen
delta_e1 = eps1_s2(1)-eps1_e(1);

%Found the elastic modulus in the plots 
%E_meas - elastic modulus from cross-head measurements
%E_corr - elastic modulus from extensometer
E_meas = 24000;
E_corr = 210000;

%Finds the corrected strain eps_ce

eps_ce = @(eps_meas, delta_eps,sigma_ce) ...
    eps_meas-delta_eps-((E_corr-E_meas)/(E_corr*E_meas))*sigma_ce;

eps_ce1= eps_ce(eps1_s2,delta_e1,sigma_ce);

%Plots the corrected stress-strain relationship of test 1:
figure();
hold on;
grid on;
plot(eps_ce1,sigma_ce, 'LineWidth',2);
plot(eps1_e,sigma_ce,'LineWidth',2);
title('Corrected stress-strain relationship, test 1');
xlabel('Engineering strain');
ylabel('Engineering stress, [MPa]');
leg3 = legend("Cross-head", "Extensometer");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task3_stress_strain_1.png');
hold off;

%% Task 4-True stress and strain
eps_log = log(1+eps_ce1);
sigma_t = sigma_ce.*(1+eps_ce1);
E = 210000;%[MPa]

%Plots the true stress-strain relationship:
figure();
hold on;
grid on;
plot(eps_log,sigma_t);
plot(eps_ce1,sigma_ce);
title('Cauchy stress-strain vs. engineering stress-strain relationship');
xlabel('Strain');
ylabel('Stress, [MPa]');
leg3 = legend("Cauchy stress", "Engineering stress");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task4_true_engineering_stress_strain.png');
hold off;

%finds the true logarithmic plastic strain
eps_log_p = eps_log - sigma_t/E;

%Plots the true stress-strain relationship:
figure();
hold on;
grid on;
plot(eps_log_p,sigma_t);
title('Cauchy stress-strain relationship');
xlabel('Plastic strain');
ylabel('Cauchy stress, [MPa]');
saveas(gcf,'task4_true_stress_strain.png');
hold off;

%% 5 - Fitted relations

%Fitted relation of the true stress
%Coefficients:
  p1 = -5.8022e+09;
  p2 = -6.0876e+09;
  p3 = 1.0754e+10;
  p4 = -5.4851e+09;
  p5 = 1.4012e+09;
  p6 = -2.0151e+08;
  p7 = 1.6619e+07;
  p8 = -7.5922e+05;
  p9 = 18527;
  p10 = 158.63;

y = @(x) p1*x^9 + p2*x^8 +...
      p3*x^7 + p4*x^6 +...
      p5*x^5 + p6*x^4 +...
      p7*x^3 + p8*x^2 +...
      p9*x + p10; 
  
eps_t = linspace(eps_log_p(1), eps_log_p(length(eps_log_p)),10);

sigma_fit = zeros(1,length(eps_t));
for i= 1:1:length(eps_t)
    sigma_fit(i) = y(eps_t(i));
end

T=table(transpose(sigma_fit),transpose(eps_t));

%% 5 - Power law and Voce law

%defining constants
sigma0_p = 286; %[MPa]
sigma0_v = 262.1;%[MPa]
K = 508.8;
n = 0.485;
Q1 = 222.7;
C1 = 10.535;

%Defining the R-function
R_p = @(p) K*p^n;

sigma_eq_p = zeros(1,length(eps_log_p));

for j = 1:1:length(R_power)
    sigma_eq_p(j)= sigma0_p + R_p(eps_log_p(j));
end


% sigma_eq_v = zeros(1,length(eps_log_p));
% for i = 1:1:length(sigma_eq_v)
%     if i ~= 1
%         sigma_eq_v(i) = sigma0_v + Q1*(1-exp(-eps_log_p(i))) + sigma_eq_v(i-1);
%     else
%         sigma_eq_v(i) = sigma0_v +Q1*(1-exp(-eps_log_p(i)));
%     end
% end

%Plots the true stress-strain relationship:
figure();
hold on;
grid on;
plot(eps_log_p,sigma_t);
plot(eps_log_p, sigma_eq_p);
%plot(eps_log_p, sigma_eq_v);
title('Cauchy stress-strain relationship');
xlabel('Plastic strain');
ylabel('Stress, [MPa]');
leg3 = legend("Cauchy stress", "Voce law","Power law");
set(leg3, 'Location', 'SouthEast');
saveas(gcf,'task5_true_stress_strain.png');
hold off;
