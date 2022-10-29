% Загрузка данных из Excel (биомасса)
% real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
% real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
% real_data(:,1) = real_data(:,1).*1000;
data = [];

for j = 3:3:44
rangeDate = strcat('C', num2str(j), ':C', num2str(j+2));
rangeX1 = strcat('D', num2str(j), ':D', num2str(j+2));
rangeX2 = strcat('F', num2str(j), ':F', num2str(j+2));
rangeStation = strcat('B', num2str(j), ':B', num2str(j+2));

Tabledata = readtable("data_groups.xlsx",'ReadVariableNames',false, 'Range', rangeDate);
rows = size(Tabledata, 1);
date = datetime(Tabledata.Var1);
diff_date = caldays(caldiff(date, 'days'));
DATE = [];
for i = 1:3:length(diff_date)
    DATE = [DATE; 0; diff_date(i); diff_date(i) + diff_date(i+1)];
end
% Загрузка данных из Excel (численность)
station = xlsread("data_groups.xlsx", rangeStation);
real_data(:,1) = xlsread("data_groups.xlsx", rangeX1);
real_data(:,2) = xlsread("data_groups.xlsx", rangeX2);
real_data(:,1) = real_data(:,1)./100000;
real_data(:,2) = real_data(:,2)./1000;

x1(1) = real_data(1, 1); %начальные данные фитопланктона
x2(1) = real_data(1, 2); %начальные данные зоопланктона

mas = [];
M1 = zeros(1, 3);
M2 = zeros(1, 3);
MIN = [10^6 10^6 10^6];

N = DATE(end); %время
h = 1; %шаг
M = 0:h:N; %сетка времени
e1 = 10^6; %критерий стабильности
flag = 0;

for alpha1 = 0.01:0.01:1
for alpha2 = 0.01:0.01:1
for beta1 = 0.01:0.01:1
for beta2 = 0.01:0.01:1
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
        %нахождение размера ошибки
%         for i=1:3
         M1 = (real_data(3, 1) - x1(DATE(3, :) + 1))^2;
         M2 = (real_data(3, 2) - x2(DATE(3, :) + 1))^2;
%         end
        %сравнение рассчитанной ошибки с минимальной
        if M1 + M2 < MIN
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
    data(end+1,:) = [station(1) mas(end,:)];
end