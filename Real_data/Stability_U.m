N = 154; %период
h = 1; %шаг

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U(1) = 0;

e1 = 10^6; %критерий устойчивости
flag = 0;
time = [0; 84; 153];
real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
real_data(:,1) = real_data(:,1)./100000;
real_data(:,2) = real_data(:,2)./1000;

x1(1) = real_data(1, 1); %начальные условия численности фитопланктона
x2(1) = real_data(1, 2); %начальные условия численности зоопланктона

alpha1(1) = 0.5; %начальные условия питания
alpha2 = 0.08;
beta1 = 0.018;
beta2 = 0.229;

mas = [];
data = [];

for k = 2:20
      MIN = 10^6;
    xc = real_data(1, 1).*k; %целевое значение
    e2 = .01 * xc; %критерий достижения цели
    mas = [];
for T1 = 0:.1:100
for T2 = 0:.1:100
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
%   период цветения воды с 1 июня (10 отчет) по 31 августа (100 отчет)
    if flag == 0 && T1 ~= 0 && T2 ~= 0
    for i=1:100
       if abs(x1(i) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && abs(x1(100) - xc) <= e2 && i < MIN
            MIN = i;
            mas(end + 1, :) = [k i T1 T2 alpha1(1) alpha2 beta1 beta2];
            break;
       end
    end
    end

    flag = 0;
end
end
     MIN = 10^6;
     MIN_mas = [];
     %находим в полученном массиве mas минимальное время
     for i = 1:size(mas,1)
           if mas(i,2) < MIN
              MIN = mas(i,2);
              MIN_mas = mas(i,:);
           end
     end
     if ~isempty(MIN_mas)
         data(end+1,:) = MIN_mas;
     end

end
%полученный на выходе массив data содержит значения коэффициентов для
%каждого значения k