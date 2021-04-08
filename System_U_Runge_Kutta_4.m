clear;

%% Классический метод Рунге-Кутта четвертого порядка

N = 100; %интервал времени
h = 1; %шаг
M = 0:h:N; %сетка времени
x1 = zeros(1, length(M)); %массив популяции жертв
x2 = zeros(1, length(M)); %массив популяции хищников
alpha1 = zeros(1, length(M)); %массив значений питания
alpha1(1) = 0.04; %начальное значение питания
x1(1) = 3.0; %начальное значение жертв
x2(1) = 4.0; %начальное значение хищников
k = 0.82; %кратность увеличения жертв
xc = k * x1(1); %целевое значение
e2 = .01 * xc; %критерий достижения цели
U(1) = 0; %управление

%коэффициенты
T1 = 1.4;
T2 = 14;
alpha2 = 0.12;
beta1 = 0.13;
beta2 = 0.05;

%функции системы ДУ
f1 = @(alpha1,x1,x2) alpha1*x1 - beta1*x1*x2;
f2 = @(x1,x2) -alpha2*x2 + beta2*x1*x2;
dfdt = @(alpha1,x1,x2) - (xc / (T2*x1*x1)) * f1(alpha1,x1,x2) + beta1 * f2(x1,x2);
fi = @(x1,x2) -( (x1 - xc)/(T2*x1) ) + beta1*x2;
psi1 = @(alpha1,x1,x2) alpha1 - fi(x1,x2);
f3 = @(alpha1,x1,x2) -(psi1(alpha1,x1,x2) / T1) + dfdt(alpha1,x1,x2);

for n=1:length(M) - 1 
    k1 = f1(alpha1(n), x1(n), x2(n));
    q1 = f2(x1(n), x2(n));
    z1 = f3(alpha1(n), x1(n), x2(n));
    
    k2 = f1(alpha1(n) + h/2*z1, x1(n) + h/2*k1, x2(n) + h/2*q1);
    q2 = f2(x1(n) + h/2*k1, x2(n) + h/2*q1);
    z2 = f3(alpha1(n) + h/2*z1, x1(n) + h/2*k1, x2(n) + h/2*q1);
    
    k3 = f1(alpha1(n) + h/2*z2, x1(n) + h/2*k2, x2(n) + h/2*q2);
    q3 = f2(x1(n) + h/2*k2, x2(n) + h/2*q2);
    z3 = f3(alpha1(n) + h/2*z2, x1(n) + h/2*k2, x2(n) + h/2*q2);
    
    k4 = f1(alpha1(n) + h*z3, x1(n) + h*k3, x2(n) + h*q3);
    q4 = f2(x1(n) + h*k3, x2(n) + h*q3);
    z4 = f3(alpha1(n) + h*z3, x1(n) + h*k3, x2(n) + h*q3);

    x1(n+1) = x1(n) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    x2(n+1) = x2(n) + h/6*(q1 + 2*q2 + 2*q3 + q4);
    U(n+1) = h/6*(z1 + 2*z2 + 2*z3 + z4);
    alpha1(n + 1) = alpha1(n) + U(n);
end

plot(M, x1, 'g','Linewidth',3);
axis([0 N -inf inf]);
hold on;
plot(M, x2, 'r', 'Linewidth',3);
hold on;
plot(M, alpha1, 'b', 'Linewidth', 3);
hold on;
plot(M, ones(length(M)).*xc, '--k','Linewidth',2);
xlabel('Время, дни');
ylabel('Популяция, ед/л');
legend('Жертва', 'Хищник', 'Питание', 'Целевое значение', 'Location','best');
ax = gca;
ax.FontSize = 20;