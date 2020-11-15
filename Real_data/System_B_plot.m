N = 154; %период
h = 1; %шаг

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U(1) = 0;

flag = 0;
real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
real_data(:,1) = real_data(:,1)./100000;
real_data(:,2) = real_data(:,2)./1000;

k = 10;
xc = real_data(1, 1).*k;
x1(1) = real_data(1, 1); %начальные условия численности фитопланктона
x2(1) = real_data(1, 2); %начальные условия численности зоопланктона
alpha1(1) = 0.5; %начальные условия питания

alpha2 = 0.0001;
beta1 = 0.01;
beta2 = 0.00033;
B = 20.0;
T1 = 3.0;
T2 = 33688.0;

for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc; 
       P = 4 * f1 * (1 / (cosh(psi)*cosh(psi)));
       dPdt = 4 * (((U(n) * x1(n) + alpha1(n) * f1 - beta1 * (f1 * x2(n) + x1(n) * f2)) * (exp(psi) + exp(-psi))^2 - 2 * f1^2 * (exp(2 * psi) - exp(-2 * psi))) / (exp(psi) + exp(-psi))^4);
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + B* P *f1;
       dfidt = ((-(T2^(-1) * dpsi0dt + df2dt) * B * P * x1(n) + B * (dPdt * x1(n) + P * f1) * (T2^(-1) * psi0 + f2)) / (B * P * x1(n))^2) + beta1 * f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -((T1)^(-1) * psi1) + dfidt;
       f3 = U(n);
       
       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
end

plot(M, x1, 'g','Linewidth',3);
axis([0 N -inf inf]);
hold on;
plot(M, x2, 'r', 'Linewidth',3);
hold on;
plot(M, alpha1, 'b', 'Linewidth', 3);
hold on;
plot(M, ones(length(M)).*xc, '--k','Linewidth',2);
hold on;
plot(M, ones(length(M)).*B, '--m','Linewidth',2);
text = strcat('Увеличение жертв в',num2str(k), ' раз');
title(text);
xlabel('Время, дни');
ylabel('Численность, экз/л');
legend('Фитопланктон (жертва)', 'Зоопланктон (хищник)', 'Питание', 'Цель', 'B');