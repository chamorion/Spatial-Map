% standard data cleaning
function [time, x, y, spk_time] = std_clean(time, x, y, spk_time)
vld = (logical(1-isnan(x))) & (logical(1-isnan(y))) & (logical(1-isnan(time)));
x = x(vld);
y = y(vld);
time = time(vld);
x_center = sum(x)/length(x);
y_center = sum(y)/length(y);
vld = logical(abs(x - x_center)<=35) & logical(abs(y - y_center)<=35);
x = x(vld);
y = y(vld);
time = time(vld);
x = x - min(x);
y = y - min(y);
vld = logical(spk_time <= max(time)) & logical(spk_time >= min(time));
spk_time = spk_time(vld);
end