function tt = read_fred_csv(csvFile)
% Read FRED csv (DATE,value). Converts '.' or missing to NaN.
T = readtable(csvFile, 'Delimiter',',', 'ReadVariableNames', true);

% DATE -> datetime
dt = datetime(T{:,1}, 'InputFormat','yyyy-MM-dd', 'TimeZone','UTC');
val = T{:,2};
if iscell(val)
    % sometimes readtable returns cellstr
    val = str2double(val);
end
val = double(val);

tt = timetable(dt, val, 'VariableNames', {'value'});
end
