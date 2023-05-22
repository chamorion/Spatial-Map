classdef Rate_Matrix
    properties
        cell_id;
        time;
        x;
        y;
        spk_time;
        itp_x;
        itp_y;
        spk_x;
        spk_y;
        rate_map;
    end
        
    methods
        function rm = Rate_Matrix(cell_id, time, x, y, spk_time)
            rm.cell_id = cell_id;
            [rm.time, rm.x, rm.y, rm.spk_time] = lib.Rate_Matrix.std_clean(time, x, y, spk_time);
            
            % interpolant function of xy coordinates vs time
            rm.itp_x=griddedInterpolant(rm.time,rm.x);
            rm.itp_y=griddedInterpolant(rm.time,rm.y);
            rm.spk_x=rm.itp_x(rm.spk_time);
            rm.spk_y=rm.itp_y(rm.spk_time);
        end
        
        % plot trajactory of the rat
        function plot_trajactory(rm)
            figure
            hold on
            title('trajactory')
            legend('show')
            plot(rm.x,rm.y,'b-','DisplayName','path of real data')
            plot(rm.spk_x,rm.spk_y,'r.','DisplayName','spikes')
        end
        
        function cal_rate_map(rm)
        end
    end
    
    methods (Static)
        % standard data cleaning
        function [time, x, y, spk_time] = std_clean(time, x, y, spk_time)
            % drop nan value
            vld = (logical(1-isnan(x))) & (logical(1-isnan(y))) & (logical(1-isnan(time)));
            x = x(vld);
            y = y(vld);
            time = time(vld);

            % centering the data and drop invalid outliers
            x_center = mean(x);
            y_center = mean(y);
            vld = logical(abs(x - x_center)<=35) & logical(abs(y - y_center)<=35);
            x = x(vld);
            y = y(vld);
            time = time(vld);
            x = x - min(x);
            y = y - min(y);

            % drop invalid time of spikes
            vld = logical(spk_time <= max(time)) & logical(spk_time >= min(time));
            spk_time = spk_time(vld);
        end
    end
    
end