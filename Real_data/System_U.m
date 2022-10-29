% clear all;
rng(1);
N = 100; %время
h = .1; %step
M = 0:h:N; %time grid

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
PSI= zeros(1, length(M));
% time = [0; 84; N];
% real_data(:,1) = xlsread("data.xlsx", 'H40:H42');
% real_data(:,2) = xlsread("data.xlsx", 'J40:J42');
% real_data(:,1) = real_data(:,1)./100000;
% real_data(:,2) = real_data(:,2)./1000;

% x1(1) = real_data(1, 1); %начальное значение жерв
% x2(1) = real_data(1, 2); %начальное значение хищников
x1(1) = 4;
x2(1) = 20;
alpha1(1) = 12.45; %начальное значение питания
% k = .17; %кратность увеличения жертв
% xc = 8; %целевое значение
% xc = k * x1(1);
U(1) = 0;

%коэффициенты
T1 = 1;
T2 = 1;
alpha2 = 0.12;
beta1 = 0.13;
beta2 = 0.05;
xc = alpha2/beta2;
e2 = .01 * xc; %критерий достижения цели

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

% tiledlayout(2,2);
% nexttile()
subplot(221);
plot(M, alpha1, 'b', 'Linewidth', 3);
hold on;
xlabel('Время, дни');
ylabel('Количество, ед/л');
legend('Питание', 'Location','best');
ax = gca;
ax.FontSize = 20;

% nexttile()
subplot(223);
plot(M, x1, 'g','Linewidth',3);
axis([0 N -inf inf]);
hold on;
% plot(M, x2, 'r', 'Linewidth',3);
% hold on;
str = strcat("\alpha_{1}(0) = ", num2str(alpha1(1)));
text(20,50,str,'FontSize',16);
plot(M, ones(length(M)).*xc, '--k','Linewidth',2);
hold on;
xlabel('Время, дни');
ylabel('Популяция, ед/л');
legend('Жертва', 'Целевое значение', 'Location','best');
ax = gca;
ax.FontSize = 20;
%% Добавление белого шума
clear all;

N = 100; %время
h = .1; %step
M = 0:h:N; %time grid
mas = [];

x = [4 20 0.04]; %начальные значения
U(1) = 0;

%коэффициенты
[T1, T2, alpha2, beta1, beta2] = deal(1, 1, 0.12, 0.13, 0.05);

% for beta2 = 0.01:0.01:1
    xc = alpha2/beta2;
for sko=0.01:0.01:10
    rng(1);
    noise = normrnd(0, sko, [1,size(M,2)]);
for n=1:(length(M) - 1)
    f1 = x(n,3)*x(n,1) - beta1*x(n,1)*x(n,2);
    f2 = -alpha2*x(n,2) + beta2*x(n,1)*x(n,2);

    dfdt = - (xc / (T2 * x(n,1) * x(n,1))) * f1 + beta1 * f2;
    fi = -( (x(n,1) - xc)/(T2*x(n,1)) ) + beta1*x(n,2);
    psi1 = x(n,3) - fi;
    U(n + 1) = -(psi1 / T1) + dfdt;
    f3 = U(n) + noise(n);

    x(n + 1,1) = x(n,1) + h*f1;
    x(n + 1,2) = x(n,2) + h*f2;
    x(n + 1,3) = x(n,3) + h*f3;
end
    x1_without_outlier = rmoutliers(x(:,1));
    e2_1 = xc - 2*std(x1_without_outlier); %левая граница
    e2_2 = xc + 2*std(x1_without_outlier); %правая граница
    flag_xc = 0; % счетчик для отслеживания достижения цели
    for j=1:length(M)
       if ~(x(j,1) >= e2_1 && x(j,1) <= e2_2)
            flag_xc = flag_xc + 1;
       end
    end
    if flag_xc / length(M) <= 0.1
        mas(end + 1,:) = [x(1,1) x(1,2) x(1,3) sko T1 T2 alpha2 beta1 beta2];
    end
end
% end

%% построение графика с шумом
rng(1);
N = 100;
h = .1;
M = 0:h:N;

i = 245;
row = num2cell(mas(i,:));
% [x(1,1), x(1,2), x(1,3), sko, T1, T2, alpha2, beta1, beta2] = deal(row{:});
[x(1,1), x(1,2), x(1,3), sko, T1, T2, alpha2, beta1, beta2] = deal(row{1:3},3.48, row{5:end});
xc = alpha2/beta2;
m_x1 = 0;
noise = normrnd(m_x1, sko, [1,size(M,2)]);

for n=1:(length(M) - 1)
    f1 = x(n,3)*x(n,1) - beta1*x(n,1)*x(n,2);
    f2 = -alpha2*x(n,2) + beta2*x(n,1)*x(n,2);

    dfdt = - (xc / (T2 * x(n,1) * x(n,1))) * f1 + beta1 * f2;
    fi = -( (x(n,1) - xc)/(T2*x(n,1)) ) + beta1*x(n,2);
    psi1 = x(n,3) - fi;
    U(n + 1) = -(psi1 / T1) + dfdt;
    f3 = U(n) + noise(n);

    x(n + 1,1) = x(n,1) + h*f1;
    x(n + 1,2) = x(n,2) + h*f2;
    x(n + 1,3) = x(n,3) + h*f3;
    
    if ~(x(n,1) >= e2_1 && x(n,1) <= e2_2)
        t = M(n);
        tx = x(n,1);
    end
end
x1_without_outlier = rmoutliers(x(:,1));
e2_1 = xc - 2*std(x1_without_outlier); %левая граница
e2_2 = xc + 2*std(x1_without_outlier); %правая граница

tiledlayout(2,2);
nexttile();
plot(M, x(:,3), 'b', 'Linewidth', 3);
xlabel('Время, дни');
ylabel('Количество, ед/л');
legend('Питание', 'Location','best');
str = strcat("\sigma = ", num2str(sko));
str1 = strcat("\mu = ", num2str(m_x1));
str = {"Параметры белого шума", str, str1};
text(5,8,str,'FontSize',16);
ax = gca;
ax.FontSize = 20;

criteria = [ones(length(M),1).*xc ones(length(M),1).*e2_1 ones(length(M),1).*e2_2];
nexttile();
hold on;
plot(M, x(:,1), 'g','Linewidth',3);
plot(M, x(:,2), 'r', 'Linewidth',3);
axis([0 N -inf inf]);
plot(M, criteria(:,1), 'k--','Linewidth',1);
plot(M, criteria(:,2),'m--','Linewidth',1);
% plot(t, tx, 'r*', 'MarkerSize', 15,'Linewidth',2);
plot(M, criteria(:,3), 'm--','Linewidth',1);
xlabel('Время, дни');
ylabel('Популяция, ед/л');
legend({'Жертва', 'Хищник', 'Целевое значение', 'Доверительный интервал'}, 'Location','best');
ax = gca;
ax.FontSize = 20;

% nexttile();
% plot(M, noise, 'black');
% axis([0 N -inf inf]);
% hold on;
% xlabel('Время, дни');
% legend('Шум', 'Location','best');
% ax = gca;
% ax.FontSize = 20;

%% генерация узкополостного случайного процесса из белого шума с фильтром
clear;
rng(1);
noise = randn(1,1000);
f = figure;
plot(noise);
title('Исходный сигнал');
ax = gca;
ax.FontSize = 20;

for i=1:1000
   noise = smooth(noise);
end

f1 = figure;
plot(noise);
title('Сигнал после 1000 прохождений фильтра smooth');
ax = gca;
ax.FontSize = 20;

f2 = figure;
histogram(noise);
ax = gca;
ax.FontSize = 20;

N = 100; %время
h = .1; %step
M = 0:h:N; %time grid

x = [4 20 0.04]; %начальные значения
U(1) = 0;

%коэффициенты
T1 = 1;
T2 = 1;
alpha2 = 0.12;
beta1 = 0.13;
beta2 = 0.05;
xc = alpha2/beta2;

for n=1:(length(M) - 1)
        f1 = x(n,3)*x(n,1) - beta1*x(n,1)*x(n,2);
        f2 = -alpha2*x(n,2) + beta2*x(n,1)*x(n,2);

        dfdt = - (xc / (T2 * x(n,1) * x(n,1))) * f1 + beta1 * f2;
        fi = -( (x(n,1) - xc)/(T2*x(n,1)) ) + beta1*x(n,2);
        psi1 = x(n,3) - fi;
        U(n + 1) = -(psi1 / T1) + dfdt;
        f3 = U(n) + noise(n);
 
        x(n + 1,1) = x(n,1) + h*f1;
        x(n + 1,2) = x(n,2) + h*f2;
        x(n + 1,3) = x(n,3) + h*f3;
end

f4 = figure;
tiledlayout(2,2);
nexttile();
plot(M, x(:,3), 'b', 'Linewidth', 3);
xlabel('Время, дни');
ylabel('Количество, ед/л');
legend('Питание', 'Location','best');
str = "1000 прохождений фильтра";
text(5,1,str,'FontSize',16);
ax = gca;
ax.FontSize = 20;

nexttile();
plot(M, x(:,1), 'g','Linewidth',3);
axis([0 N -inf inf]);
hold on;
plot(M, x(:,2), 'r', 'Linewidth',3);
hold on;
plot(M, ones(length(M)).*xc, '--k','Linewidth',2);
xlabel('Время, дни');
ylabel('Популяция, ед/л');
legend('Жертва', 'Хищник', 'Целевое значение', 'Location','best');
ax = gca;
ax.FontSize = 20;