%Problem set 6- task 2 - Newton Raphson method
% initializing constants
Q = 100; 
C = 10;
H = 50; 
F = 200;
tol = 0.001;

% creating functions that handles F and dF/dx
F_int = @(x) Q.*(1-exp(-C.*x)) + H.*x;
dF_int = @(x) Q.*C.*exp(-C.*x) + H;
x = 0.2205;
x_n = x;
n = 0;

%initializing a matrix with the iteration values
M(1,1) = n;
M(1,2) = x_n;
M(1,3) = -(F_int(x_n)-F)./(dF_int(x_n));
M(1,4) = abs(F_int(x_n)-F);

% loop iterates until one is below the tolerance
while abs(F_int(x_n)-F) > tol
 %The residual G is given by F_int(x_n)-F
 x_n = x_n - (F_int(x_n)-F)./(dF_int(x_n));
 n = n+1;
 M(n+1,1) = n;
 M(n+1,2) = x_n;
 M(n+1,3) = -(F_int(x_n)-F)./(dF_int(x_n));
 M(n+1,4) = abs(F_int(x_n)-F);
 M(n,5) = (M(n+1,2)-M(n,2))./M(n,2);
end