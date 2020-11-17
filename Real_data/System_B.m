N = 154; %период
h = 1; %шаг

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U(1) = 0;

e1 = 10^3; %критерий устойчивости
flag = 0;
real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
real_data(:,1) = real_data(:,1)./100000;
real_data(:,2) = real_data(:,2)./1000;

k = 25;
B = 20;
xc = real_data(1, 1).*k;
e2 = .01 * xc; %критерий достижения цели
x1(1) = real_data(1, 1); %начальные условия численности фитопланктона
x2(1) = real_data(1, 2); %начальные условия численности зоопланктона
x1(1) = x1(1) - xc / 2;
x2(1) = x2(1) - B / 2;

mas = [];
mas_stab = [];
data = [];
max_stab = 0;
flag_find = 0;
T2_min = 0;
T2_max = 1000000;

% alpha1(1) = 0.12;
% alpha2 = 0.06;
% beta1 = 0.002;
% beta2 = 0.082;

for a1 = 0.1:0.1:1
    alpha1(1) = a1;
for alpha2 = 0.01:0.01:1
for beta1 = 0.0001:0.0001:0.1
for beta2 = 0:0.00005:0.001
% for B = 100:10:200
for T1 = 0:1:100000
for T2 = T2_min:T2_max
for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc; 
       P = 4 * f1 * (1 / (cosh(psi)*cosh(psi)));
       dPdt = 4 * (((U(n) * x1(n) + alpha1(n) * f1 - beta1 * (f1 * x2(n) + x1(n) * f2)) * (exp(psi) + exp(-psi))^2 - 2 * f1^2 * (exp(2 * psi) - exp(-2 * psi))) / (exp(psi) + exp(-psi))^4);
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + B* P *f1;
       dfidt = ((-(T2^(-1) * dpsi0dt + df2dt) * B * P * x1(n) + B * (dPdt * x1(n) + P * f1) * (T2^(-1) * psi0 + f2)) / (B * P * x1(n))^2) + beta1 * f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -((T1)^(-1) * psi1) + dfidt;
       f3 = U(n);

       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
       
        if abs(x1(n + 1)) >= 2*xc
           flag = 1;
           if n > max_stab
              max_stab = n; 
              disp([T1 T2]);
           end
            break
         end
        if abs(x2(n + 1)) > B
           flag = 1;
           if n > max_stab
              max_stab = n; 
              disp([T1 T2]);
           end
           break
        end
         if abs(alpha1(n + 1)) >= e1
            flag = 1;
            if n > max_stab
              max_stab = n; 
              disp([T1 T2]);
           end
            break
         end
end
    if flag == 0 && T1 ~= 0 && T2 ~= 0 && beta2 ~= 0
        mas_stab(end + 1, :) = [k B T1 T2 alpha1(1) alpha2 beta1 beta2];
        flag_find = 1;
        disp([k B T1 T2 alpha1(1) alpha2 beta1 beta2]);
%     for i=1:100
%        if abs(x1(i) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && abs(x1(100) - xc) <= e2
%             mas(end + 1, :) = [k i T1 T2 alpha1(1) alpha2 beta1 beta2];
%             disp([k B i T1 T2 alpha1(1) alpha2 beta1 beta2]);
%             break;
%        end
%     end
    end
    flag = 0;
    if flag_find == 1
       break; 
    end
end
    if flag_find == 1
        flag_find = 0;
%         T2_min = T2_min + 1000;
        T2_max = T2_max + 5000;
       break; 
    end
end
end
end
end
end
% end
% end