%проверка Т1 и Т2
load T.mat;
h = 0.3;
alpha2 = 0.2;
beta1 = 0.8;
beta2 = 0.1;
M = 0:h:N; %сетка времени
x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
x1(1) = 1; %начальные условия численности жертв
x2(1) = 1; %начальные условия численности хищников
alpha1(1) = 1; %начальные условия питания
U(1) = 0;
for j = 1:1:size(T, 1)
    for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

       dfdt = - (xc / (T(j, 2) * x1(n) * x1(n))) * f1 + beta1 * f2;
       fi = -( (x1(n) - xc)/(T(j, 2)*x1(n)) ) + beta1*x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -(psi1 / T(j, 1)) + dfdt;

       f3 = U(n);

       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
    end
    if((j > 0) && (j <= 15))
       subplot(321);
    end
    if((j > 15) && (j <= 30))
       subplot(322);
    end
    if((j > 30) && (j <= 45))
       subplot(323);
    end
    if((j > 45) && (j <= 60))
       subplot(324);
    end
    if((j > 60) && (j <= 71))
       subplot(325);
    end
    plot(M, x1, 'g', M, x2, 'r', M, alpha1, 'b');
       xlabel('Время');
       ylabel('Популяция');
       legend('Жертвы', 'Хищники', 'Питание');
       hold on;
end