load('2.mat'); %load file with data
new_H_T1_T2 = zeros(1, 3); %new massiv with correct value
h_T2 = zeros(1, 2); %massiv with h and T2 (min)
for i = 1:size(H_T1_T2,1)
    if H_T1_T2(i, 2) > H_T1_T2(i, 1)
        new_H_T1_T2(end + 1, :) = H_T1_T2(i, :);
%         if H_T1_T2(i, 2) - 0.01 >= H_T1_T2(i, 1)
%            h_T2(end + 1, :) = [H_T1_T2(i, 1) H_T1_T2(i, 2)];
%         end
    end
end
for i = 1:size(H_T1_T2,1)-1
    if H_T1_T2(i+1, 1) > H_T1_T2(i, 1)
       h_T2(end + 1, :) = [H_T1_T2(i+1, 1) H_T1_T2(i+1, 2)];
    end
end
% plot3(new_H_T1_T2(:,1), new_H_T1_T2(:,2), new_H_T1_T2(:,3),'o');
% xlabel('h');
% ylabel('T2');
% zlabel('T1');
% subplot(221);
% plot(new_H_T1_T2(2:100,2), new_H_T1_T2(2:100,3));
% text = strcat('h = ',num2str(new_H_T1_T2(2,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(222);
% plot(new_H_T1_T2(392:486,2), new_H_T1_T2(392:486,3));
% text = strcat('h = ',num2str(new_H_T1_T2(392,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(223);
% plot(new_H_T1_T2(857:946,2), new_H_T1_T2(857:946,3));
% text = strcat('h = ',num2str(new_H_T1_T2(857,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(224);
% plot(new_H_T1_T2(1297:1381,2), new_H_T1_T2(1297:1381,3));
% text = strcat('h = ',num2str(new_H_T1_T2(1297,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% plot3(new_H_T1_T2(:,1), new_H_T1_T2(:,2), new_H_T1_T2(:,3));
% xlabel('h');
% ylabel('T2');
% zlabel('T1');
[c error] = polyfit(h_T2(:,1),h_T2(:,2),2);
[Y delta] = polyval(c,h_T2(:,1),error);
% plot(h_T2(:,1),Y,'g',h_T2(:,1),Y+2*delta,'r-',h_T2(:,1),Y-2*delta,'r-')
% title('График h, T2 и теоретическая функция');
% hold on;
% plot(h_T2(:,1),h_T2(:,2),'o');
% plot(new_H_T1_T2(2:100,2),polyval(c,new_H_T1_T2(2:100,2)));
% text = strcat('h = ',num2str(new_H_T1_T2(2,1)));
% title(text);
% xlabel('h');
% ylabel('T2');