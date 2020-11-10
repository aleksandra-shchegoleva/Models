N = 154; %period of the time
h = 1; %step

M = 0:h:N; %time grid

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));

time = [0; 85; 154];
real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
real_data(:,1) = real_data(:,1).*1000;
real_data(:,2) = real_data(:,2);
x1(1) = real_data(1, 1); %initial conditions of phytoplankton
x2(1) = real_data(1, 2); %initial conditions of zooplankton
xc = real_data(3, 1); %goal value
M1 = zeros(1, 3);
M2 = zeros(1, 3);

alpha1 = 0.03;
alpha2 = 0.87;
beta1 = 0.001;
beta2 = 0.01;

for n=1:(length(M) - 1)
      f1 = alpha1*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

       x1(n + 1) = x1(n) + h*f1;
        x2(n + 1) = x2(n) + h*f2;
end

plot(M, x1, 'g', 'Linewidth',4);
hold on;
plot(M, x2, 'r', 'Linewidth',4), 
hold on;
% plot(time(3,1), xc, 'g*', 'MarkerSize', 20);
axis([0 N -inf inf]);
hold on;
plot(time, real_data(:,1), 'g--o', 'Linewidth',2);
hold on;
plot(time, real_data(:,2), 'r--o', 'Linewidth',2);
xlabel('Время, дни');
ylabel('Численность, экз/л');
legend('Фитопланктон (жертва)', 'Зоопланктон (хищник)','Реальные данные (фитопланктон)','Реальные данные (зоопланктон)');
D_x1 = var(x1);
D_x2 = var(x2);
S_x1 = std(x1);
S_x2 = std(x2);
stat = [D_x1 D_x2 S_x1 S_x2];
num = 1;
            mas_max = [];
            MAX = 0;
            while(num ~= 0)
                MAX = 0;
                num = 0;
                for n=2:(length(M) - 1)
                    if(x1(n) > x1(n-1) && x1(n) > x1(n+1) && x1(n) > MAX && x1(n) ~= 10^6)
                        MAX = x1(n);
                        num = n;
                    end
                end
                if(num ~= 0)
                    x1(num) = 10^6;
                    mas_max(:, end + 1) = MAX;
                end
            end