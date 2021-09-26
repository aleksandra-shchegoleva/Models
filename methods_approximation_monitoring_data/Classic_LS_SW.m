%% Нахождение коэффициентов alpha1, alpha2, beta1, beta2 для дискретной системы "хищник-жертва" 
%% по данным мониторинга с помощью классического метода наименьших квадратов

clc;
clear all;

% данные загружаются в виде таблицы
data = readtable('SW_2019_data\2019_B1_phytoplankton_zooplankton_aphanocapsa_cyclopoida.csv');
x_data = [data.ValuePhytoplankton data.ValueZooplankton];
date = data.Sampling_date; % дата представлена в формате YYYY-MM-DD

% расчет разности между датами
t_diff = days(diff(date));
t = 0;
% расчет номера дня с начала наблюдений для каждого измерения
for i=1:length(t_diff)
    t(end+1) = t(i) + t_diff(i);
end

h = 1; % шаг дискретизации
% все переменные рациональные числа, сделать их положительными - заменить
% real на positive (возможно решение не будет найдено)
syms alpha1 alpha2 beta1 beta2 real;
[dQda1, dQda2, dQdb1, dQdb2] = deal(0,0,0,0); % сумма производных по каждому параметру
for i=1:size(x_data,1)-1
    dQda1 = dQda1 + ((x_data(i,1) + h*(alpha1*x_data(i,1) - beta1*x_data(i,1)*x_data(i,2))) - x_data(i+1,1))*(h*x_data(i,1)); 
    dQda2 = dQda2 + ((x_data(i,2) + h*(-alpha2*x_data(i,2) + beta2*x_data(i,1)*x_data(i,2))) - x_data(i+1,2))*(-h*x_data(i,2));
    dQdb1 = dQdb1 + ((x_data(i,1) + h*(alpha1*x_data(i,1) - beta1*x_data(i,1)*x_data(i,2))) - x_data(i+1,1))*(-h*x_data(i,1)*x_data(i,2));
    dQdb2 = dQdb2 + ((x_data(i,2) + h*(-alpha2*x_data(i,2) + beta2*x_data(i,1)*x_data(i,2))) - x_data(i+1,2))*(h*x_data(i,1)*x_data(i,2));
end

% требуется найти решение, при котором сумма квадратов отклонений будет
% равна 0
dQ = [dQda1 == 0 dQda2 == 0 dQdb1 == 0 dQdb2 == 0];

[a1, a2, b1, b2] = solve(dQ,[alpha1 alpha2 beta1 beta2]); % решение системы dQ
[a1, a2, b1, b2] = deal(double(a1),double(a2),double(b1),double(b2)); % приведение к десятиной дроби

disp(strcat('Станция: ',string(data.Reported_station_name(1))));
disp(strcat('Фитопланктон: ', string(data.Scientific_name(1)), ' Зоопланктон: ', string(data.Scientific_name_1(1))));
disp(strcat('Классический МНК: a1=', string(a1),' b1=', string(b1),' a2=', string(a2), ' b2=', string(b2)));

%% построение графика

x = zeros(1,2);
[x(1,1), x(1,2)] = deal(x_data(1,1), x_data(1,2));
M = 0:h:t(end);

for n=1:(length(M) - 1)
    f1 = a1*x(n,1) - b1*x(n,1)*x(n,2);
    f2 = -a2*x(n,2) + b2*x(n,1)*x(n,2);

    x(n + 1,1) = x(n,1) + h*f1;
    x(n + 1,2) = x(n,2) + h*f2;
end

figure;
% окно 2x2
tiledlayout(2,2);
nexttile
plot(date,x_data,'Linewidth',3);
xlim([date(1) date(end)]);
title('Данные мониторинга');
legend({'Жертвы (данные мониторинга)','Хищники (данные мониторинга)'});
ax = gca;
ax.FontSize = 20;

nexttile
hold on;
plot(M,x(:,1) ,'g--','Linewidth',3);
plot(M,x(:,2) ,'r--','Linewidth',3);
title('Классический МНК');
legend({'Жертвы (функция)','Хищники (функция)'});
ax = gca;
ax.FontSize = 20;

% наложение рассчитанных значений на данные мониторинга
nexttile([1 2])
hold on;
plot(date,x_data,'Linewidth',3);
plot(M,x(:,1) ,'g--','Linewidth',3);
plot(M,x(:,2) ,'r--','Linewidth',3);
xlim([date(1) date(end)]);
title('Наложение рассчитанных значений на данные мониторинга');
legend({'Жертвы (данные мониторинга)','Хищники (данные мониторинга)','Жертвы (функция)','Хищники (функция)'});
ax = gca;
ax.FontSize = 20;