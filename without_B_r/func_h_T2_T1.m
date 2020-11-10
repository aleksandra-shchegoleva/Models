load('3-1.mat'); %load file with data
new_H_T1_T2 = zeros(1, 3); %new massiv with correct value
h_T2 = zeros(1, 2); %massiv with h and T2 (min)
h_T1 = zeros(1, 2); %massiv with h and T1 (min)
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
for i = 1:size(H_T1_T2,1)-1
    if H_T1_T2(i+1, 1) > H_T1_T2(i, 1)
       h_T1(end + 1, :) = [H_T1_T2(i+1, 1) H_T1_T2(i+1, 3)];
    end
end
% plot3(new_H_T1_T2(:,1), new_H_T1_T2(:,2), new_H_T1_T2(:,3), 'o');
% xlabel('h');
% ylabel('T2');
% zlabel('T1');
% subplot(221);
% plot(new_H_T1_T2(2:99,2), new_H_T1_T2(2:99,3));
% text = strcat('h = ',num2str(new_H_T1_T2(2,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(222);
% plot(new_H_T1_T2(376:461,2), new_H_T1_T2(376:461,3));
% text = strcat('h = ',num2str(new_H_T1_T2(376,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(223);
% plot(new_H_T1_T2(772:837,2), new_H_T1_T2(772:837,3));
% text = strcat('h = ',num2str(new_H_T1_T2(772,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% subplot(224);
% plot(new_H_T1_T2(1054:1091,2), new_H_T1_T2(1054:1091,3));
% text = strcat('h = ',num2str(new_H_T1_T2(1054,1)));
% title(text);
% xlabel('T2');
% ylabel('T1');
% plot3(new_H_T1_T2(:,1), new_H_T1_T2(:,2), new_H_T1_T2(:,3));
% xlabel('h');
% ylabel('T2');
% zlabel('T1');
[c error] = polyfit(h_T1(:,1),h_T1(:,2),2);
[Y delta] = polyval(c,h_T1(:,1),error);
% plot(h_T1(:,1),Y,'g',h_T1(:,1),Y+2*delta,'r-',h_T1(:,1),Y-2*delta,'r-')
% title('График h, T1 и теоретическая функция');
% hold on;
% plot(h_T1(:,1),h_T1(:,2),'o');
% plot(new_H_T1_T2(2:100,2),polyval(c,new_H_T1_T2(2:100,2)));
% text = strcat('h = ',num2str(new_H_T1_T2(2,1)));
% title(text);
% xlabel('h');
% ylabel('T2');