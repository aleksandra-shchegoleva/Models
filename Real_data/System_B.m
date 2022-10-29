N = 154; %�����
h = 1; %���

M = 0:h:N; %����� �������

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U = zeros(1, length(M));
U(1) = 0;

e1 = 10^2; %�������� ������������
flag = 0;
real_data(:,1) = xlsread("data.xlsx", 'H3:H5');
real_data(:,2) = xlsread("data.xlsx", 'J3:J5');
real_data(:,1) = real_data(:,1)./100000;
real_data(:,2) = real_data(:,2)./1000;

% k = 1/10;
% B = 30;
% xc = 10;
% r = 0.01;
x1(1) = real_data(1, 1); %��������� ������� ����������� �������������
x2(1) = real_data(1, 2); %��������� ������� ����������� ������������
% x1(1) = 10;
% x2(1) = 1;
% x1(1) = x1(1) - xc / 2;
% x2(1) = x2(1) - B / 2;

mas = [];
mas_stab = [];
data = [];
max_stab = 0;
flag_find = 0;
T2_min = 1;
T2_max = 1000;
R = [];

alpha1(1) = 0.13; %��������� ������� �������
alpha2 = 0.2;
beta1 = 0.1;
beta2 = 0.3;
% T1 = 2.537;
% T2 = 11.951;
% r = 0.001;

% for a1 = 0.1:0.1:10
%     alpha1(1) = a1;
% for alpha2 = 0.001:0.001:1
% for beta1 = 0.0001:0.0001:0.01
% for beta2 = 0.0001:0.0003:0.01
for k = 0.1:0.1:10
    xc = x1(1) * k;
    e2 = .01 * xc; %�������� ���������� ����
for B = 1:1:5
% for r = 0.01:0.01:1
for T1 = 0:0.1:10
for T2 = 0:1:100
for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
       
       psi = x1(n) - xc; 
       P = 1 / (cosh(psi))^2;
       dPdt = -2*(cosh(psi))^-3 * sinh(psi) * f1;
%        P = 4 * f1 * (1 / (exp(psi) + exp(-psi))^2);
%        dPdt = 4 * (((U(n) * x1(n) + alpha1(n) * f1 - beta1 * (f1 * x2(n) + x1(n) * f2)) * (exp(psi) + exp(-psi))^2 - 2 * f1^2 * (exp(2 * psi) - exp(-2 * psi))) / (exp(psi) + exp(-psi))^4);
       df2dt = -alpha2 * f2 + beta2 * (f1 * x2(n) + x1(n) * f2);
       psi0 = x2(n) + B*tanh(psi);
       dpsi0dt = f2 + B* P *f1;
       dfidt = -(((T2^-1 * dpsi0dt + df2dt)*B*P*x1(n) - B*(dPdt*x1(n) + P*f1))*(T2^-1*psi0 + f2))/(B*P*x1(n))^2 + beta2*f2;
       fi = - (T2 ^ (-1) * psi0 + f2) / (B * P * x1(n)) + beta1 * x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -((T1)^(-1) * psi1) + dfidt;
       r = abs(psi*P/xc);
       f3 = U(n);
       U(n+1) = r*U(n+1);

       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
       
        if abs(x1(n + 1)) >= e1 || isnan(x1(n+1))
           flag = 1;
           if n > max_stab
              max_stab = n; 
%               disp([k B i T1 T2 alpha1(1) alpha2 beta1 beta2]);
           end
            break
         end
        if abs(x2(n + 1)) > e1 || isnan(x2(n+1))
           flag = 1;
           if n > max_stab
              max_stab = n; 
%               disp([k B i T1 T2 alpha1(1) alpha2 beta1 beta2]);
           end
           break
        end
        if abs(alpha1(n + 1)) > e1 || isnan(alpha1(n+1))
           flag = 1;
           if n > max_stab
              max_stab = n; 
%               disp([T1 T2]);
           end
           break;
        end
        
        if isnan(U(n+1)) || abs(U(n+1)) > e1
           flag = 1;
           break;
        end
end
    if flag == 0 && T1 ~= 0 && T2 ~= 0
        mas_stab(end + 1, :) = [xc B T1 T2 alpha1(1) alpha2 beta1 beta2];
%         disp([k B T1 T2 alpha1(1) alpha2 beta1 beta2]);
        for i=1:length(M)
           if abs(x1(i) - xc) <= e2 && abs(mean(x1(i:end)) - xc) <= e2 && abs(x1(end) - xc) <= e2 && abs(x1(end-10) - xc) <= e2%&& abs(x1(i + 1) - xc) <= e2 && abs(x1(end-40) - xc) <= e2
                mas(end + 1, :) = [k B T1 T2 alpha1(1) alpha2 beta1 beta2];
                flag_find = 1;
%                 disp([k B i T1 T2 alpha1(1) alpha2 beta1 beta2]);
                break;
           end
        end
    end

    flag = 0;
    if flag_find == 1
        break; 
    end
end
    if flag_find == 1
        flag_find = 0;
        break; 
    end
end
% if flag_find == 1
%         flag_find = 0;
% % %         T2_min = T2_min + 1000;
% % %         T2_max = T2_max + 5000;
%        break; 
%     end
end
end
% end
% end
% end
% end
% end