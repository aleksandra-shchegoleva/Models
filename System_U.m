N = 200; %максимум сетки
h = 0.02; %шаг

alpha2 = 0.05; %коэффициент смертности хищников
beta1 = 0.4; %коэффициент смертности жертв
beta2 = 0.01; %коэффициент рождаемости хищников
 
xc = 5.0; %целевое значение
 
T2 = 5.0;
T1 = 0.5;

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
x1(1) = 20.0; %начальные условия численности жертв
x2(1) = 5.0; %начальные условия численности хищников
alpha1(1) = 0; %начальные условия питания
U(1) = 0;

for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

        dfdt = - (xc / (T2 * x1(n) * x1(n))) * f1 + beta1 * f2;
        fi = -( (x1(n) - xc)/(T2*x1(n)) ) + beta1*x2(n);
        psi1 = alpha1(n) - fi;
        U(n + 1) = -(psi1 / T1) + dfdt;
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
xlabel('Время, дни');
ylabel('Биомасса, мг/л');
legend('Фитопланктон (жертва)', 'Зоопланктон (хищник)', 'Питание', 'Цель');