%Backward euler method with work hardening
%The algorithm is defined as: 
%y_(n+1) = y_n + f(t_(n+1),y_(n+1))*delta_t_(n+1)

%constants
E = 210000; %Elastic modulus [Mpa]
h_R = E/9; %hardening modulus [Mpa]
sigma0 = E/1000; %yield stress [Mpa]
alpha = 0.4;
epsilon0 = 1/1000;
epsilon1 = (1+alpha)*epsilon0;
t0 = 1;

%creating vectors
t = linspace(0,t0,1000); %time [s]
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
       epsilon(n) = epsilon1.*(t(n)./t-4);
   end
end

