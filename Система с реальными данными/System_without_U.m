N = 154; %period of the time
h = .1; %step

M = 0:h:N; %time grid

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));

time = [0; 85; 154];
real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
real_data(:,1) = real_data(:,1).*1000;
x1(1) = real_data(1, 1); %initial conditions of phytoplankton
x2(1) = real_data(1, 2); %initial conditions of zooplankton
xc = real_data(3, 1); %goal value
M1 = zeros(1, 3);
M2 = zeros(1, 3);

alpha1 = 0.0001;
alpha2 = 0.0008;
beta1 = 0.034;
beta2 = 0.002;

for n=1:(length(M) - 1)
      f1 = alpha1*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

       x1(n + 1) = x1(n) + h*f1;
        x2(n + 1) = x2(n) + h*f2;
end

plot(M, x1, 'g', M, x2, 'r', time(3,1), xc, 'g*', 'MarkerSize', 20);
axis([0 N -inf inf]);
hold on;
plot(time, real_data(:,1), 'g--o', time, real_data(:,2), 'r--o');
xlabel('Время, дни');
ylabel('Биомасса, мг/куб м');
legend('Фитопланктон (жертва)', 'Зоопланктон (хищник)');
D_x1 = var(x1);
D_x2 = var(x2);
S_x1 = std(x1);
S_x2 = std(x2);
stat = [D_x1 D_x2 S_x1 S_x2];
for i=1:3
   M1(1, i) = (real_data(i, 1) - x1(time(i, :) + 1))^2;
   M2(1, i) = (real_data(i, 2) - x2(time(i, :) + 1))^2;
end