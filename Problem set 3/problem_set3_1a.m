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
        for n = 1:1:N
            R_temp = R(p(n));
            sigma_vec(x,n) = sigma(m_temp,dp_temp,R_temp);
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

    