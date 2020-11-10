N = 154; %время
h = 1; %шаг

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));

e1 = 10^6; %критерий стабильности
flag = 0;
time = [0; 85; 154];

choose = input("Выбери формат исходных данных: биомасса (1), популяция (2)\n");
if choose == 1
    real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
    real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
    real_data(:,1) = real_data(:,1).*1000;
    real_data(:,2) = real_data(:,2);
    labely = "Биомасса, мг/л";
else
    real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
    real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
    real_data(:,1) = real_data(:,1)./100000;
    real_data(:,2) = real_data(:,2)./1000;
    labely = "Популяция, экз/м3";
end
k = input("УВеличение жертв ");

x1(1) = real_data(1, 1); %начальное значение фитопланктона
x2(1) = real_data(1, 2); %начальное значение зоопланктона

mas = [];

M1 = zeros(1, 3);
M2 = zeros(1, 3);
MIN = [10^6 10^6 10^6];

for alpha1 = 0.01:0.01:0.5
for alpha2 = 0.01:0.01:0.5
for beta1 = 0.001:0.001:0.5
for beta2 = 0.001:0.001:0.5
    
for n=1:(length(M) - 1)
      f1 = alpha1*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
 
       x1(n + 1) = x1(n) + h*f1;
        x2(n + 1) = x2(n) + h*f2;

        if x1(n + 1) >= e1 || x1(n + 1) < 0
           flag = 1;
           break
        end
        if x2(n + 1) >= e1 || x2(n + 1) < 0
           flag = 1;
           break
        end
end
    if flag == 0
        for i=1:3
         M1(1, i) = (real_data(i, 1) - x1(time(i, :) + 1))^2;
         M2(1, i) = (real_data(i, 2) - x2(time(i, :) + 1))^2;
        end
        if sum(M1) + sum(M2) < sum(MIN)
            MIN = M1 + M2;
            mas(end + 1, :) = [alpha1 alpha2 beta1 beta2];
        end
    end
    M1 = [];
    M2 = [];
    flag = 0;
end
end
end
end

alpha1 = zeros(1, length(M));
U(1) = 0;

xc = x1(1) * k; %целевое значение
e2 = 0.05 * xc;
alpha1(1) = mas(1); %начальное значение питания
%массив mas хранит полученные коэффициенты
alpha2 = mas(2);
beta1 = mas(3);
beta2 = mas(4);
MIN = 100;
data = [];

for j = size(mas,1):-1:floor(size(mas,1)/2)
    alpha1(1) = mas(j,1);
    alpha2 = mas(j,2);
    beta1 = mas(j,3);
    beta2 = mas(j,4);
for T1 = 0:100
for T2 = 0:100
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
        if x1(n + 1) >= e1 || x1(n + 1) < 0
           flag = 1;
           break
        end
        if x2(n + 1) >= e1 || x2(n + 1) < 0
           flag = 1;
           break
        end
        if alpha1(n + 1) >= e1  || alpha1(n + 1) < 0
           flag = 1;
           break
        end
end

    if flag == 0 && T1 ~= 0 && T2 ~= 0
    for i=1:100
       if abs(x1(i) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && abs(x1(100) - xc) <= e2&& i < MIN
            MIN = i;
            data = [i T1 T2 alpha1(1) alpha2 beta1 beta2];
            break;
       end
    end
    end

    flag = 0;
end
end
end


T1 = data(2);
T2 = data(3);
alpha1(1) = data(4);
alpha2 = data(5);
beta1 = data(6);
beta2 = data(7);

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
text = strcat('Увеличение by',num2str(k), ' times');
title(text);
xlabel('Time, days');
ylabel(labely);
legend('Phytoplankton (prey)', 'Zooplankton (predator)', 'Food', 'Goal');