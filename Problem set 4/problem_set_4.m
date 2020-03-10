%Problem set 4: Tensile test

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
%%



