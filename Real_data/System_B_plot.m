N = 154; %������
h = 1; %���

M = 0:h:N; %����� �������

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
alpha1 = zeros(1, length(M));
U(1) = 0;

flag = 0;
time = [0; 84; 153];
real_data(:,1) = xlsread("data.xlsx", 'I3:I5');
real_data(:,2) = xlsread("data.xlsx", 'K3:K5');
real_data(:,1) = real_data(:,1).*1000;
real_data(:,2) = real_data(:,2);

k = 1.5;
xc = real_data(1, 1).*k;
x1(1) = real_data(1, 1); %��������� ������� ����������� �������������
x2(1) = real_data(1, 2); %��������� ������� ����������� ������������
alpha1(1) = 0.5; %��������� ������� �������

alpha2 = 0.01;
beta1 = 0.0021;
beta2 = 0.0011;
B = 500.0;
T1 = 1.9;
T2 = 1925000;

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
end

plot(M, x1, 'g','Linewidth',3);
axis([0 N -inf inf]);
hold on;
plot(M, x2, 'r', 'Linewidth',3);
hold on;
plot(M, alpha1, 'b', 'Linewidth', 3);
hold on;
plot(M, ones(length(M)).*xc, '--k','Linewidth',2);
hold on;
plot(M, ones(length(M)).*B, '--m','Linewidth',2);
text = strcat('���������� ����� �',num2str(k), ' ���');
title(text);
xlabel('�����, ���');
ylabel('��������, ��/�');
legend('������������ (������)', '����������� (������)', '�������', '����', 'B');