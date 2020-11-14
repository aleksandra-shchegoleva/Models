Tprop = readtable("PROP_stability_u.csv");
prop = table2array(Tprop);
N = prop(1); %максимум сетки
h = prop(2); %шаг

T = readtable("DATA_stability_u.csv");
data = table2array(T);
k = input('Введите номер строчки из массива с коэффициентами дл¤ построени¤ графика');
alpha2 = data(k, 5); %коэффициент смертности хищников
beta1 = data(k, 6); %коэффициент смертности жертв
beta2 = data(k, 7); %коэффициент рождаемости хищников
 
xc = prop(5); %целевое значение
 
T1 = data(k,2);
T2 = data(k,3);

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
x1(1) = prop(3); %начальные условия численности жертв
x2(1) = prop(4); %начальные условия численности хищников
alpha1(1) = data(k,4); %начальные условия питани¤
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