function x = transform_yoy_log(level_series)
% YoY log growth: log(x_t) - log(x_{t-12})
x = log(level_series) - log(lagmatrix(level_series,12));
end
