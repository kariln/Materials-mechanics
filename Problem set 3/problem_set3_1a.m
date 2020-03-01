%Problem set 3: viscoplasticity and large deformations
%task 1a

%defining variables
m = linspace(0.02,0.1,5); %rate sensitivities
K = 500; %strength coefficient[Mpa] 
dp0 = 1/1000; %reference strain rate [s^-1]
n = 0.3; %hardening exponent
dp = [1/1000,1,1000]; %plastic strain rates

%defining vectors
N = 1000; %number of sampling points
p = linspace(0,0.05,N);
sigma_vec = zeros(length(m),N);

%defines functions
R = @ (p)...
    K*p^n;
sigma = @ (m,dp,R) ...
    (dp/dp0)^m*R;

%updating values
for y = 1:1:length(dp)
    dp_temp = dp(y);
    figure();
    grid on;
    for x = 1:1:length(m)
        m_temp = m(x);
        for z = 1:1:N
            R_temp = R(p(z));
            sigma_vec(x,z) = sigma(m_temp,dp_temp,R_temp);
        end
        txt = ['m = ',num2str(m_temp)];
        plot(p,sigma_vec(x,:),'DisplayName',txt);
        hold on;
        title("Plastic strain rate = "+ dp_temp);
        xlabel('Equivalent plastic strain, p [mm/mm]');
        ylabel('Equivalent plastic stress, sigma [MPa]');
        legend show;
        file_name = dp_temp +".png";
        saveas(gcf,file_name)
    end
end

%%
%task 1d
%plotting p(t) with several values of m

N = 1000;
t = linspace(0,100,N);
p = zeros(length(m),N); 
load = 200; %[Mpa]

%creating functions
p_func = @(t,m)... %the equivalent plastic strain as a function of t
    ((m+n)/m*dp0*(load/K)^(1/m)*t)^(m/(n+m));

%calculating the values with the function
figure();
hold on;
grid on;
xlabel('Time, t [s]');
ylabel('Equivalent plastic strain, p')
title('Equivalent plastic strain as a function of time');
for i = 1:1:length(m)
    m_tmp = m(i);
    for j = 1:1:N
        p(i,j) = p_func(t(j),m_tmp); 
    end
    text = ['m = ',num2str(m_tmp)];
    plot(t,p(i,:),'DisplayName',text);
    hold on;
end
legend show;
saveas(gcf,'eq_pl_str.png');
hold off;
