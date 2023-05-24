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
            [rm.time, rm.x, rm.y, rm.spk_time] = lib.Rate_Matrix.std_preprocess(time, x, y, spk_time);
            
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
            rm.rate_map=flipud(rot90(hm_spk./hm_time));
        end              
        
        % spatial smoothing based on Gaussian kernel
        function rm = spatial_smooth(rm, bin_num)
            delta_x = (max(rm.x) - min(rm.x))/bin_num;
            delta_y = (max(rm.y) - min(rm.y))/bin_num;
            delta_t = mean(diff(rm.time));
            sigma = min([delta_x, delta_y])/3;
            Gauss_K = @(v1, v2) exp(- (v1.^2 + v2.^2)/2*sigma^2);
            
            t = min(rm.time) : delta_t : max(rm.time);
            p_x = @(t) rm.itp_x(t);
            p_y = @(t) rm.itp_y(t);

            rm.rate_map = zeros(bin_num, bin_num);
            i_x = 1;
            for s_x = min(rm.x) + delta_x/2 : delta_x : max(rm.x) - delta_x/2
                i_y = 1;
                for s_y = min(rm.y) + delta_y/2 : delta_y : max(rm.y) - delta_y/2
                    rm.rate_map(i_x, i_y) = sum(Gauss_K(rm.spk_x - s_x, rm.spk_y - s_y))/trapz(Gauss_K(p_x(t) - s_x, p_y(t) - s_y));
                    i_y = i_y + 1;
                end
                i_x = i_x + 1;
            end
            rm.rate_map = flipud(rot90(rm.rate_map));
            hm_xy=histcounts2(rm.x,rm.y,'NumBins',[bin_num, bin_num]);
            rm.rate_map(flipud(rot90(hm_xy)) == 0)=NaN;
        end
        
        % plot trajactory of the cell
        function plot_trajactory(rm)
            figure;
            hold on;
            title(['trajactory of cell ',num2str(rm.cell_id)]);
            legend('show');
            plot(rm.x,rm.y,'b-','DisplayName','trajactory');
            plot(rm.spk_x,rm.spk_y,'r.','DisplayName','spikes');
        end
        
        % plot place field of the cell
        function plot_rate_map(rm)
            figure;
            hold on;
            title(['place field of cell ',num2str(rm.cell_id)]);
            hm=imagesc(rm.rate_map);
            set(gca,'YDir','normal');
            colormap('turbo');
            set(hm,'AlphaData',~isnan(hm.CData));
            colorbar;
        end
    end
    
    methods (Static)
        % standard data preprocessing: data cleaning and centering
        function [time, x, y, spk_time] = std_preprocess(time, x, y, spk_time)
            % drop nan value
            vld = (logical(1-isnan(x))) & (logical(1-isnan(y)));
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
            
            % x,y centering
            x = x - min(x);
            y = y - min(y);
            
            % drop invalid time of spikes
            vld = logical(spk_time <= max(time)) & logical(spk_time >= min(time));
            spk_time = spk_time(vld);
            
            % time centering
            spk_time = spk_time - min(time);
            time = time - min(time);
        end
        
    end
    
end