N = 500; %time
h = 0.01;
M = 0:h:N; %time grid
 
x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U = zeros(1, length(M));

alpha2 = 0.38; %predator attrition rate
beta1 = 0.01; %prey attrition rate
beta2 = 0.005; %predator fertility rate
 
x1(1) = 10; %initial conditions of population preys
x2(1) = 2; %initial conditions of population predators
alpha1(1) = 0; %initial conditions of food
U(1) = 0;
xc = 100.0;
B = 15; %population of predators
 
T1 = 1.0;
T2 = 1.0;

for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc; 
       P = 4 / (exp(psi) + exp(-psi))^2;
        r = abs(psi * P  / xc);
       dPdt = (8*(exp(-psi)*f1 - exp(psi)*f1))/(exp(-psi) + exp(psi))^3;
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + B* P *f1;
       dfidt = ((-(T2^(-1) * dpsi0dt + df2dt) * B * P * x1(n) + B * (dPdt * x1(n) + P * f1) * (T2^(-1) * psi0 + f2)) / (B * P * x1(n))^2) + beta1 * f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -((T1)^(-1) * psi1) + dfidt;
       f3 = r*U(n);
       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
end
  plot(M, x1, 'g', M, x2, 'r', M, alpha1, 'b');
  xlabel('Время');
  ylabel('Популяция');
  legend('Жертвы', 'Хищники', 'Питание');