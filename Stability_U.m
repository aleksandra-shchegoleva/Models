N = input("Введите время существования системы\n");
h = input("Введите шаг\n");
M = 0:h:N; %сетка времени
x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));

xc = input("Введите целевое значение\n");
x1(1) = input("Начальные условия численности жертв\n");
x2(1) = input("Начальные условия численности хищников\n");
alpha1(1) = input("Начальные условия питания\n");

e1 = 10^6; %критерий для устройчивости системы
e2 = .05*xc; %критерий для достижения цели
flag = 0;
MIN = 10^6;
data = [];
U(1) = 0;

for alpha2 = 0.1:0.1:1
for beta1 = 0.001:0.001:0.1
for beta2 = 0.01:0.01:0.5
for T1 = 1:2:50
for T2 = 70:2:100
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
    if abs(x1(n + 1)) >= e1
        flag = 1;
        break;
    end
    if abs(x2(n + 1)) >= e1
        flag = 1;
        break;
    end
    if abs(alpha1(n + 1)) >= e1
        flag = 1;
        break;
    end
end
    if flag == 0 
       for i=1:length(M) - 1
           if abs(x1(i) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && abs(x1(length(M)) - xc) <= e2 && i < MIN
                MIN = i;
                data(end + 1,:) = [i T1 T2 alpha1(1) alpha2 beta1 beta2];
                break;
           end
        end 
    end
end
end
end
end
end
varNamesP = {'Time', 'h', 'x1', 'x2', 'xc'};
varNamesC = {'TimeTarget', 'T1', 'T2', 'alpha1', 'alpha2', 'beta1', 'beta2'};
properties = [N h x1(1) x2(1) xc];
Tprop = array2table(properties, 'VariableNames', varNamesP);
T = array2table(data, 'VariableNames', varNamesC);
writetable(T, "DATA_stability_u.csv");
writetable(Tprop, "PROP_stability_u.csv");