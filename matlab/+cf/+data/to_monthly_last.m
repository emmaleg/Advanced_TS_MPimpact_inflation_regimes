function mt = to_monthly_last(tt)
% Aggregate (possibly daily) timetable to monthly by taking last non-missing observation.

% Remove NaNs first (keep last non-NaN)
tt2 = tt;
tt2 = tt2(~isnan(tt2.value),:);

% If already monthly (one obs per month), retime still safe.
mt = retime(tt2, 'monthly', 'lastvalue');
end
