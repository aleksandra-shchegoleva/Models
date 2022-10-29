function y=plot_system_full(N, x1Init, x2Init, xc, date)

h = 0.01;
M = 0:h:N; %сетка времени

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
x_data = [x1Init, x2Init];

% начальные значения коэффициентов
coeff = [0.01 0.01 0.01 0.01];

% нелинейный МНК (Квазиньютоновский метод)
F = @(coeff) func(coeff,x_data(1,:),date); % расчет значений в точках
FsumSquares = @(coeff) sum(sum((F(coeff) - x_data).^2))/size(x_data,1); % критерий МНК
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[quasiCoeff,ressquared,eflag,outputu] = fminunc(FsumSquares,coeff,opts);
[alpha1, beta1, alpha2, beta2] = deal(quasiCoeff(1), quasiCoeff(2), quasiCoeff(3), quasiCoeff(4));
clc;

e1 = 10^6; %критерий стабильности
flag = 0;

x1(1) = x1Init(1); %начальное значение фитопланктона
x2(1) = x2Init(1); %начальное значение зоопланктона

mas = [];

alpha1 = zeros(1, length(M));
alpha1(1) = quasiCoeff(1);
u(1) = 0;

disp(quasiCoeff)
if xc == -1
    xc = alpha2/beta2; %целевое значение
end
e2 = 0.05 * xc;
MIN = 100;
data = [];

% for T1 = 0:.5:50
% for T2 = 0:.5:50
% for n=1:(length(M) - 1)
%       f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
%        f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
% 
%         dfdt = - (xc / (T2 * x1(n) * x1(n))) * f1 + beta1 * f2;
%         fi = -( (x1(n) - xc)/(T2*x1(n)) ) + beta1*x2(n);
%         psi1 = alpha1(n) - fi;
%         U(n + 1) = -(psi1 / T1) + dfdt;
%         f3 = U(n);
%  
%        x1(n + 1) = x1(n) + h*f1;
%         x2(n + 1) = x2(n) + h*f2;
%         alpha1(n + 1) = alpha1(n) + h*f3;
%         if x1(n + 1) >= e1 || x1(n + 1) < 0
%            flag = 1;
%            break
%         end
%         if x2(n + 1) >= e1 || x2(n + 1) < 0
%            flag = 1;
%            break
%         end
%         if alpha1(n + 1) >= e1  || alpha1(n + 1) < 0
%            flag = 1;
%            break
%         end
% end
% 
%     if flag == 0 && T1 ~= 0 && T2 ~= 0
%     for i=1:100
%        if abs(x1(i) - xc) <= e2 && abs(x1(i + 1) - xc) <= e2 && abs(x1(100) - xc) <= e2&& i < MIN
%             MIN = i;
%             data = [i T1 T2 alpha1(1) alpha2 beta1 beta2];
%             break;
%        end
%     end
%     end
% 
%     flag = 0;
% end
% end

% if ~isempty(data)
%     T1 = data(2);
%     T2 = data(3);
%     alpha1(1) = data(4);
%     alpha2 = data(5);
%     beta1 = data(6);
%     beta2 = data(7);
T1 = 5; T2 = 5;
    for n=1:(length(M) - 1)
          f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
           f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

            dfdt = - (xc / (T2 * x1(n) * x1(n))) * f1 + beta1 * f2;
            fi = -( (x1(n) - xc)/(T2*x1(n)) ) + beta1*x2(n);
            psi1 = alpha1(n) - fi;
            u(n + 1) = -(psi1 / T1) + dfdt;
            f3 = u(n);

           x1(n + 1) = x1(n) + h*f1;
            x2(n + 1) = x2(n) + h*f2;
            alpha1(n + 1) = alpha1(n) + h*f3;
    end
    y = [x1; x2; alpha1];
% else
%    y = 'Система не имеет решений'; 
% end
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