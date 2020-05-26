e1 = 1000; %criterion for system stability
xc = 10.0;
e2 = .2*xc; %criteria for a goal
flag = 0;
h = 0.1;
N = 10; %time

U(1) = 0;
B = xc; %population of predators
alpha2 = 0.9; %predator attrition rate
beta1 = 0.1; %prey attrition rate
beta2 = 0.1; %predator fertility rate
coeff = zeros(1, 2);

% for h = 0.1:0.1:1.0
    M = 0:h:N; %time grid
    x1 = zeros(1, length(M));
    x2 = zeros(1, length(M));
    alpha1 = zeros(1, length(M));
    x1(1) = 5; %initial conditions of population preys
    x2(1) = 2; %initial conditions of population predators
    alpha1(1) = 2; %initial conditions of food
% for alpha2 = 0.1:0.1:1.5
% for beta1 = 0.1:0.1:1.5
% for beta2 = 0.1:0.1:1.5
for T1 = -50:1:0
for T2 = -50.0:1:50
for n = 1:(length(M) - 1)
    f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc;
       P = 4 / (exp(psi) + exp(-psi))^2;
       dPdt = -2*f1*P*(exp(2 * psi) - exp(-2 * psi));
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + P*f1;
       dfidt = ((-(T2^(-1) * dpsi0dt + df2dt) * B * P * x1(n) + B * (dPdt * x1(n) + P * f1) * (T2^(-1) * psi0 + f2)) / (B * P * x1(n))^2) + beta1 * f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -(T1^(-1) * psi1) + dfidt;

       f3 = U(n);

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
       coeff(end + 1, :) = [T1 T2];
    end
    flag = 0;
end
end
% end
% end
% end
% end

% save('coeff.mat', 'coeff');