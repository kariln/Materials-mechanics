%Backward euler method with work hardening
%The algorithm is defined as: 
%y_(n+1) = y_n + f(t_(n+1),y_(n+1))*delta_t_(n+1)

%constants
E = 210000; %Elastic modulus [Mpa]
h_R = E/9; %hardening modulus [Mpa]
sigma0 = E/1000; %yield stress [Mpa]
alpha = 0.4;
epsilon0 = 1/1000;
epsilon1 = (1+alpha).*epsilon0;
t0 = 1;

%creating vectors
t = linspace(0,4*t0,1000); %time [s]
N = length(t);
epsilon = zeros(1,N); %strain history, epsilon(t)
p = epsilon; %accumulated plastic strain
R = epsilon; %hardening variable, R = R(p)
epsilon_e = epsilon; %elastic strain
epsilon_p = epsilon; %plastic strain
e_p0 = 0; %initial plastic strain
p_0 = 0; %initial accumulated plastic strain

%creating yield function using function_handle:
%f = abs(sigma_tr) - (sigma_0 + R_tr)
%where sigma_tr and R_tr is stress and hardening in the trial state
f = @(sigma_tr,sigma_0,R_tr)...
    abs(sigma_tr)-(sigma_0+R_tr);

%creating increment for p, dp:
dp = @(sigma_tr, sigma_0,R_n,E,h_R)...
    (abs(sigma_tr)-(sigma_0+R_n))/(E+h_R);


%defines the strain history, epsilon(t), as it was defined in task 1b
for n = 1:1:N
   if t(n) <= t0
       epsilon(n) = t(n).*epsilon1./t0;
   elseif t(n) > t0 && t(n) <= 3*t0
       epsilon(n) = epsilon1.*(2-t(n)./t0);
   elseif t(n) > 3*t0
       epsilon(n) = epsilon1.*(t(n)./t0-4);
   end
end

%updating variables
for n = 1:1:N-1
    %creating the trial state
    sigma_tr = sigma(n) + E.*(epsilon(n+1)-epsilon(n));
    R_tr = R(n);
    
    %checks if the trial state fulfils the yield conditions at t_(n+1)
    if f(sigma_tr,sigma0,R_tr) <= 0
        %fulfils yield condition
        sigma(n+1) = sigma_tr;
        R(n+1) = R_tr;
        epsilon_p(n+1) = epsilon_p(n);
        epsilon_e(n+1) = epsilon(n+1)-epsilon_p(n+1);
        p(n+1) = p(n);
        
    else
        %does not fulfil the yield condition
        delta_p = dp(sigma_tr,sigma0,R(n),E,h_R);
        eps_p = sign(sigma_tr).*delta_p;
        p(n+1) = p(n) + delta_p;
        epsilon_p(n+1) = epsilon_p(n) + eps_p;
        epsilon_e(n+1) = epsilon(n+1) - epsilon_p(n+1);
        sigma(n+1) = sigma_tr - E.*(eps_p);
        R(n+1) = R(n) + h_R.*delta_p;
    end
end

figure();
P1 = plot(epsilon,sigma);
grid on;
set(P1,'linewidth',1.2)%changes the thickness of the graph line
fs = 13;
ylabel('Stress [MPa]','Fontsize',fs) %ylabel with customized fontsize
xlabel('Strain [mm/mm]','Fontsize',fs);
title('Backward-Euler algorithm for linear isotropic hardening')
saveas(gcf,'Backward-Euler.png')
