function y=parse_date(FILE)

Time = table();

%определение времени
Tabledata = readtable(FILE,'ReadVariableNames',false);
rows = size(Tabledata, 1);
date = datetime(Tabledata.Var3(3:rows,:));
diff_date = caldays(caldiff(date, 'days'));
DATE = [];
for i = 1:3:length(diff_date)
    DATE = [DATE; 0; diff_date(i); diff_date(i) + diff_date(i+1)];
end
DATE = table(DATE);
Time = [Time DATE];

%номера станций
Time = [Time, Tabledata.Var2(3:rows,:)];

y = Time;