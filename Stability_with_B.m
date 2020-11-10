flag = 0;
h = 0.01;
N = 1; %time
M = 0:h:N;

% alpha2 = 0.3; %predator attrition rate
% beta1 = 0.4; %prey attrition rate
% beta2 = 0.03; %predator fertility rate
 
% r = 0.01;
T2 = 5.0;
T1 = 5.0;
xc = 30;
B = xc; %population of predators
e1 = 10^3;
e2 = 0.95 * xc;

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
U = zeros(1, length(M));
alpha1 = zeros(1, length(M));
x1(1) = 5; %initial conditions of population preys
x2(1) = 1; %initial conditions of population predators
alpha1(1) = 0; %initial conditions of food
U(1) = 0;
mas = zeros(1, 3);

for beta2 = 0.001:0.001:1
for beta1 = 0.001:0.001:1
for alpha2 = 0.001:0.001:1
for n = 1:(length(M) - 1)
    f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc; 
       P = 4 / (exp(psi) + exp(-psi))^2;
       r = abs(psi * P / xc);
       dPdt = (8*(exp(-psi)*f1 - exp(psi)*f1))/(exp(-psi) + exp(psi))^3;
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + B* P *f1;
       dfidt = ((-(T2^(-1) * dpsi0dt + df2dt) * B * P * x1(n) + B * (dPdt * x1(n) + P * f1) * (T2^(-1) * psi0 + f2)) / (B * P * x1(n))^2) + beta1 * f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
        U(n + 1) = -((T1)^(-1) * psi1) + dfidt;
       f3 = r * U(n);
       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
    if (abs(x1(n)) >= e1)
        flag = 1;
        break;
    end
    if (abs(x2(n)) >= B*10) 
        flag = 1;
        break;
    end
    if (abs(alpha1(n)) >= e1)
        flag = 1;
        break;
    end
end
    if flag == 0
        for i = 1:length(M) - 1*(1/h)
            if abs(mean(x1(i : floor(i + 1*(1/h)))) - xc) <= e2 
                 mas(end + 1, :) = [alpha2 beta1 beta2];
                 break;
            end
        end
    end
    flag = 0;
 end
end
end
% end
% end
