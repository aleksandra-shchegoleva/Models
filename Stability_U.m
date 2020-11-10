e1 = 10^6; %критерий стабильности
xc = 10; %целевое значение
e2 = .05*xc; %критерий достижения цели
flag = 0;

N = 500; %время

U(1) = 0;
h = 0.3;
M = 0:h:N; %временная сетка
mas = [];

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
x1(1) = 5; %начальное значение жертв
x2(1) = 2; %начальное значение хищников
alpha1(1) = 0; %начальное значение питания

for alpha2 = 0.1:0.1:2.0
for beta1 = 0.1:0.1:2.0
for beta2 = 0.1:0.1:2.0
for T1 = 0.1:0.1:0.5
for T2 = 100.0:1.0:120 
    MIN = 10^6;
for n = 1:(length(M) - 1)
    f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
    f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
    dfdt = - (xc / (T2 * (x1(n) * x1(n)))) * f1 + beta1 * f2;
    fi = -( (x1(n) - xc)/(T2*x1(n)) ) + beta1*x2(n);
    psi1 = alpha1(n) - fi;
    U(n + 1) = -(psi1 / T1) + dfdt;
    f3 = U(n);

    x1(n + 1) = x1(n) + h*f1;
    x2(n + 1) = x2(n) + h*f2;
    alpha1(n + 1) = alpha1(n) + h*f3;
    
    if x1(n + 1)) >= e1 || x1(n + 1) < 0
        flag = 1;
        break;
    end
    if x2(n + 1) >= e1 || x2(n + 1) < 0
        flag = 1;
        break;
    end
    if alpha1(n + 1) >= e1 || alpha1(n + 1) < 0
        flag = 1;
        break;
    end
end
    if flag == 0 
        for i=1:100
        if abs(mean(x1(i:100)) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && i < MIN
            MIN = i;
            mas(end + 1, :) = [i T1 T2 alpha1(1) alpha2 beta1 beta2];
            break;
       end
        end
    end
    flag = 0;
end
end
end
end
end