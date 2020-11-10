e1 = 25; %критерий для устройчивости системы
xc = 10;
e2 = .2*xc; %критерий для достижения цели
flag = 0;

N = 50; %время

U(1) = 0;

fid = fopen('stability.txt', 'w'); %открытие файла для записи данных
fid1 = fopen('target.txt', 'w');

for alpha2 = 0.01:0.01:5
   for beta1 = 0.01:0.01:5
      for beta2 = 0.01:0.01:5
            for T1 = 0.1:.1:100 %величина, влияющая на длительность переходного процесса
                for T2 = 0.1:0.1:100
                    for h = 0.01:0.01:1
                         M = 0:h:N; %сетка времени

                        x1 = zeros(1, length(M));
                        x2 = zeros(1, length(M));
                        x1(1) = 1; %начальные условия численности жертв
                        x2(1) = 1; %начальные условия численности хищников
                        alpha1(1) = 1; %начальные условия питания

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
                       fprintf(fid, 'a2 = %g, beta1 = %g, beta2 = %g, T1 = %g, T2 = %g, xc = %g, h = %g\n', alpha1, beta1, beta2, T1, T2, xc, h);
                       xs = mean(x1(end - (length(M) - 20):end-2)); %среднее значение
                       if abs(xs - xc) <= e2
                            fprintf(fid1, 'T1 = %g, T2 = %g, xc = %g, h = %g\n', T1, T2, xc, h);
                       end
                 end
           flag = 0;
                    end
end
end
      end
   end
end



fclose( fid );
fclose( fid1 );