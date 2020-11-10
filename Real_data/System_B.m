N = 154; %период
h = 1; %шаг

M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U(1) = 0;

e1 = 10^3; %критерий устойчивости
flag = 0;
real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
real_data(:,1) = real_data(:,1).*1000;
real_data(:,2) = real_data(:,2);

x2(1) = real_data(1, 2); %начальные условия численности зоопланктона
alpha1(1) = 0.5; %начальные условия питания
k = 2;
xc = real_data(1, 1).*k;
e2 = .01 * xc; %критерий достижения цели
x1(1) = real_data(1, 1); %начальные условия численности фитопланктона

mas = [];
mas_stab = [];
data = [];

for alpha2 = 0.001:0.001:1
for beta1 = 0.00001:0.0001:1
for beta2 = 0.00001:0.0001:1
for B = 100:10:200
for T1 = .01:0.01:2
for T2 = 10^4:10000:10^7
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
%        data(end + 1, :) = [f1 f2 psi P dPdt dfidt fi psi1 f3];

       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
       
        if abs(x1(n + 1)) >= e1 || x1(n + 1) < 0
           flag = 1;
            break
         end
        if abs(x2(n + 1)) >= e1 || x2(n + 1) < 0
           flag = 1;
           break
        end
         if abs(alpha1(n + 1)) >= e1 || alpha1(n + 1) < 0
            flag = 1;
            break
         end
end
    if flag == 0 && T1 ~= 0 && T2 ~= 0
        mas_stab(end + 1, :) = [k B T1 T2 alpha1(1) alpha2 beta1 beta2];
%     for i=1:100
%        if abs(mean(x1(i:100)) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2
%             mas(end + 1, :) = [k i T1 T2 alpha1(1) alpha2 beta1 beta2];
%             disp([k B i T1 T2 alpha1(1) alpha2 beta1 beta2]);
%             break;
%        end
%     end
    end
    
    flag = 0;
end
end
end
end
end
end