N = 154; %врем€
h = 1; %шаг

M = 0:h:N; %сетка времени

mas = [];
x1 = zeros(1, length(M));
x2 = zeros(1, length(M));

e1 = 10^6; %stability criterion
flag = 0;
time = [0; 85; N];
% %загрузка данных из Excel файла (попул€ци€)
% real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
% real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
% real_data(:,1) = real_data(:,1)./100000;
% real_data(:,2) = real_data(:,2)./1000;
% %загрузка данных из Excel файла (биомасса)
% real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
% real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
% real_data(:,1) = real_data(:,1).*1000;
% real_data(:,2) = real_data(:,2);

x1(1) = real_data(1, 1); %initial conditions of phytoplankton
x2(1) = real_data(1, 2); %initial conditions of zooplankton
xc = real_data(3, 1); %goal value

M1 = zeros(1, 3);
M2 = zeros(1, 3);
MIN = [10^6 10^6 10^6];

for alpha1 = 0.01:0.01:0.5
for alpha2 = 0.01:0.01:1
for beta1 = 0.001:0.001:1
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
        %расчет размера ошибки
        for i=1:3
         M1(1, i) = (real_data(i, 1) - x1(time(i, :) + 1))^2;
         M2(1, i) = (real_data(i, 2) - x2(time(i, :) + 1))^2;
        end
        %проверка размера ошибки
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
