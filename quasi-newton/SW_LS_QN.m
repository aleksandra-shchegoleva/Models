% данные загружаются в виде таблицы
function [a1, b1, a2, b2] = SW_LS_QN()
data = readtable('SW_2019_data/2019_B1_phytoplankton_zooplankton_aphanocapsa_keratella_quadrata.csv');
x_data = [data.ValuePhytoplankton data.ValueZooplankton];
date = data.Sampling_date; % дата представлена в формате YYYY-MM-DD

% расчет разности между датами
t_diff = days(diff(date));
t = 0;
% расчет номера дня с начала наблюдений для каждого измерения
for i=1:length(t_diff)
    t(end+1) = t(i) + t_diff(i);
end

% начальные значения коэффициентов
coeff = [0.01 0.01 0.01 0.01];

% нелинейный МНК (Квазиньютоновский метод)
F = @(coeff) func(coeff,x_data(1,:),t); % расчет значений в точках
FsumSquares = @(coeff) sum(sum((F(coeff) - x_data).^2)); % критерий МНК
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[quasiCoeff,ressquared,eflag,outputu] = fminunc(FsumSquares,coeff,opts);
[a1, b1, a2, b2] = deal(quasiCoeff(1), quasiCoeff(2), quasiCoeff(3), quasiCoeff(4));

clc;
end

function y = func(coeff,x0,t)
    h = 1;
    x = x0;
    j = 2;
    y = x0;
    for i=1:t(end)-1
       x(i+1,1) = x(i,1) + h*x(i,1)*(coeff(1) - coeff(2)*x(i,2));
       x(i+1,2) = x(i,2) + h*x(i,2)*(-coeff(3) + coeff(4)*x(i,1));
       if i+1 == t(j)
            y(end+1, :) = x(i+1,:);
            j = j + 1;
       end
    end
end