% plot trajactory of the rat
function plot_trajactory(time, x, y, spk_time)
[time, x, y, spk_time] = std_clean(time, x, y, spk_time);
figure
hold on
title('trajactory')
legend('show')
itp_x=griddedInterpolant(time,x); % interpolant function of x coordinates vs. time coordinates
itp_y=griddedInterpolant(time,y);
spk_x=itp_x(spk_time);
spk_y=itp_y(spk_time);
plot(x,y,'b-','DisplayName','path of real data')
plot(spk_x,spk_y,'r.','DisplayName','spikes')
end