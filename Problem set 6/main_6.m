%Problem set 6
clear all;
close all;

%% Task 1b
%defining constants
E = 210000; % Young's modulus
epsilon_0 = 1.0; % Reference strain
fs = 13;

% Tolerance criterion for the palstic iterations
tol = 0.000001;

% Material constants in Voce law from Problem Set 4
sigma_0 = 262.08;
Q1 = 222.705; 
C1 = 10.535;

% Time-variables
t_0 = 1;
t_max = 1*t_0;
dt = 0.01;
t = 0:dt:t_max;

len_t = length(t);

% Strain history as defined in the problem set
epsilon = epsilon_0.*sin(pi.*t./(2.*t_0));

%Defining the yield criterion for elastic-plastic material with isotropic
%hardening
f = @(sigma,sigma_0, R) abs(sigma)-(sigma_0+R);

%Defining the plastic modulus h_r, found in 1a
h_R = @(C_R,Q_R,p) C_R*Q_R*exp(-C_R*p);

%Defining vectors of zeros
sigma = zeros(1,len_t);
R = zeros(1,len_t);
p = zeros(1,len_t);
epsilon_e = zeros(1,len_t);
epsilon_p = zeros(1,len_t);

%Defining the trial state by a loop
for n = 1:1:len_t-1
    delta_eps = epsilon(1,n+1)-epsilon(1,n);
    sigma_tr = sigma(1,n) + E*delta_eps; %elastic trial state based on the variables at t_n
    R_tr = R(1,n); %not changing during the time step as it is assumed to be elastic
    if f(sigma_0,sigma_tr,R_tr) <= 0 %if yes, the step is elastic
        sigma(1,n+1) = sigma_tr;
        R(1,n+1) = R_tr;
        p(1,n+1) = p(1,n);
        epsilon_p(1,n+1) = epsilon_p(1,n);
    else % f>0, the step is plastic, and therefore inadmissable
       %has to re-establish consistency, f=0
       %introduces temporary variables
       yield_function = f(sigma_tr,sigma_0,R_tr);
       hardening = R_tr;
       plastic = p(1,n);
       sig = sigma_tr;
       
       while yield_function > abs(tol)
        delta_p = abs(sigma_tr)-(sigma_0+R(1,n))/(E+h_R(C1,Q1,plastic));
        delta_eps_p = delta_p*sgn(sigma_tr);
        plastic = p(1,m) + delta_p;
        delta_hard = h_R(C1,Q1,plastic)*delta_p;
        hardening = R(1,n) + delta_hard;
        sig = sigma_tr - E*delta_eps_p; 
        yield_function = sgn(sigma_tr)*(sigma_tr-E*delta_eps_p)-(sigma_0+R(1,n)+h_R(C1,Q1,plastic));
       end
       
       sigma(1,n+1) = sig;
       R(1,n+1) = hardening;
       p(1,n+1) = plastic;
       epsilon_p(1,n+1) = epsilon(1,n)+delta_eps_p;
    end
    epsilon_e(1,n+1) = epsilon(1,n+1)-epsilon_p(1,n+1);
end

figure();
plot(epsilon,sigma);
xlabel('Strain');
ylabel('Stress');
title('Stress-strain curve');
axis([-0.01 0.05 0 360])
    

        

