function t =parse_date(FILE)

Time = table();

%определение времени
Tabledata = readtable(FILE,'ReadVariableNames',false);
rows = size(Tabledata, 1);
date = Tabledata.Var3;

% расчет разности между датами
t_diff = days(diff(date));
t = 0;
% расчет номера дня с начала наблюдений для каждого измерения
for i=1:length(t_diff)
    if i == 3
       t(end+1) = 0;
    else
        t(end+1) = t(i) + t_diff(i);
    end
end
