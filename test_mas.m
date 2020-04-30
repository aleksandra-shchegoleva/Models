%check T1 and T2
load T.mat;
N = 50;
xc = 10;
h = 0.3;
alpha2 = 0.2;
beta1 = 0.8;
beta2 = 0.1;
M = 0:h:N; %time grid
x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
x1(1) = 1; %initial conditions of population preys
x2(1) = 1; %initial conditions of population predators
alpha1(1) = 1; %initial conditions of food
U(1) = 0;
for j = 1:1:size(T, 1)
    for n=1:(length(M) - 1)
      f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
       f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);

       dfdt = - (xc / (T(j, 2) * x1(n) * x1(n))) * f1 + beta1 * f2;
       fi = -( (x1(n) - xc)/(T(j, 2)*x1(n)) ) + beta1*x2(n);
       psi1 = alpha1(n) - fi;
       U(n + 1) = -(psi1 / T(j, 1)) + dfdt;

       f3 = U(n);

       x1(n + 1) = x1(n) + h*f1;
       x2(n + 1) = x2(n) + h*f2;
       alpha1(n + 1) = alpha1(n) + h*f3;
    end
    if((j >= 16) && (j <= 23))
       subplot(221);
       plot(x1, x2, 'g');
       xlabel('Preys');
       ylabel('Predators');
       hold on;
       subplot(222);
       plot(x1, alpha1,'g');
       xlabel('Preys');
       ylabel('Food');
       hold on;
       subplot(223);
       plot(x2, alpha1, 'r');
       xlabel('Predators');
       ylabel('Food');
       hold on;
       subplot(224);
       plot3(x1, x2, alpha1);
       xlabel('Preys');
       ylabel('Predators');
       zlabel('Food');
       legend('Preys', 'Predators', 'Food');
       hold on;
    end
%     if((j >= 24) && (j <= 26))
%        subplot(132);
%     end
%     if((j >= 27) && (j <= 30))
%        subplot(133);
%      end
%     if((j > 45) && (j <= 60))
%        subplot(324);
%     end
%     if((j > 60) && (j <= 71))
%        subplot(325);
%     end
    
end