classdef Rate_Matrix < handle
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
        
        function rm = cal_rate_map(rm, bin_num)
            delta_t = mean(diff(rm.time));
            [hm_xy,x_edge,y_edge]=histcounts2(rm.x,rm.y,'NumBins',[bin_num, bin_num]);
            hm_xy(hm_xy==0)=NaN;
            hm_time=hm_xy*delta_t;
            hm_spk=histcounts2(rm.spk_x,rm.spk_y,'XBinEdges',x_edge,'YBinEdges',y_edge);
            rm.rate_map=hm_spk./hm_time;
        end
        
        % plot trajactory of the cell
        function plot_trajactory(rm)
            disp(rm.rate_map)
            figure;
            hold on;
            title('trajactory');
            legend('show');
            plot(rm.x,rm.y,'b-','DisplayName','path of real data');
            plot(rm.spk_x,rm.spk_y,'r.','DisplayName','spikes');
        end
        
        % plot place field of the cell
        function plot_rate_map(rm)
            figure;
            hold on;
            title('place field');
            hm=imagesc(flipud(rot90(rm.rate_map)));
            set(gca,'YDir','normal');
            colormap('turbo');
            set(hm,'AlphaData',~isnan(hm.CData));
            colorbar;
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