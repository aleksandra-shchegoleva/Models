function stations =parse_stations(FILE)

Time = table();

%определение времени
Tabledata = readtable(FILE,'ReadVariableNames',false);
stations = Tabledata.Var2;
stations = string(stations');
