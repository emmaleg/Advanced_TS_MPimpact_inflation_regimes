function tt = merge_series(tt_list, names)
% Outer-join multiple timetables on time.
% tt_list: cell array of timetables with variable 'value'
% names: variable names for each series

tt = tt_list{1};
tt.Properties.VariableNames = names(1);

for i=2:numel(tt_list)
    tmp = tt_list{i};
    tmp.Properties.VariableNames = names(i);
    tt = synchronize(tt, tmp, 'union', 'mean'); %#ok<*AGROW>
end
end
