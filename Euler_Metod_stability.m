e1 = 25; %criterion for system stability
xc = 10;
e2 = .2*xc; %criteria for a goal
flag = 0;

N = 50; %time

U(1) = 0;
h = 0.3;
alpha2 = 0.2;
beta1 = 0.8;
beta2 = 0.1;
M = 0:h:N; %time grid
T = zeros(1, 2); %table with coefficients T

x1 = zeros(1, length(M));
x2 = zeros(1, length(M));
x1(1) = 1; %initial conditions of population preys
x2(1) = 1; %initial conditions of population predators
alpha1(1) = 1; %initial conditions of food

fid = fopen('stability.txt', 'w'); %open file for data recording
fid1 = fopen('target.txt', 'w');

% for alpha2 = 0.1:0.001:0.3
    for T1 = 0.1:0.1:1
      for T2 = 80.0:1.0:110
%            for T2 = 80:1.0:110     
            for n = 1:(length(M) - 1)
                   f1 = alpha1(n)*x1(n) - beta1*x1(n)*x2(n);
                   f2 = -alpha2*x2(n) + beta2*x1(n)*x2(n);
                   dfdt = - (xc / (T2 * (x1(n) * x1(n)))) * f1 + beta1 * f2;
                   fi = -( (x1(n) - xc)/(T2*x1(n)) ) + beta1*x2(n);
                   psi1 = alpha1(n) - fi;
                   U(n + 1) = -(psi1 / T1) + dfdt;

                   f3 = U(n);

                   x1(n + 1) = x1(n) + h*f1;
                   x2(n + 1) = x2(n) + h*f2;
                   alpha1(n + 1) = alpha1(n) + h*f3;
                   if (x1(n) >= e1)
                         flag = 1;
                         break;
                   end
                   if (x2(n) >= e1) 
                         flag = 1;
                         break;
                   end
                   if (alpha1(n) >= e1)
                         flag = 1;
                         break;
                   end
            end
            if flag == 0 
                fprintf(fid, 'a2 = %g, beta1 = %g, beta2 = %g, T1 = %g, T2 = %g\n', alpha2, beta1, beta2, T1, T2);
                T(end + 1, :) = [T1 T2];
                for i = 1:length(x1)-25
                   if(abs(x1(i) - xc) <= e2) 
                       fprintf(fid1, 'a2 = %g, beta1 = %g, beta2 = %g, T1 = %g, T2 = %g\n', alpha2, beta1, beta2, T1, T2);
                       break;
                   end
                end
            end
       flag = 0;
              end
            
%            end
%            end
 end
save('T.mat', 'T');


fclose( fid );
fclose( fid1 );