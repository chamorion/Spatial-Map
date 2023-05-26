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
        bin_num;
        sp_info;
        coherence;
        sparsity;
        visits;
        unsmoothed;
        rate_map;
        auto_corr_m;
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
        
        % spatial-smoothed rate map based on Gaussian kernel
        function rm = gauss_rm(rm, bin_num)
            delta_x = (max(rm.x) - min(rm.x))/bin_num;
            delta_y = (max(rm.y) - min(rm.y))/bin_num;
            delta_t = mean(diff(rm.time));
            
            % 2D histogram of visit times in each bin
            [rm.visits, x_edge,y_edge] = histcounts2(rm.x, rm.y, 'NumBins', [bin_num, bin_num]);
            rm.visits = flipud(rot90(rm.visits));
            
            % unsmoothed rate map
            occp_time = rm.visits*delta_t;
            spk_in_bin = flipud(rot90(histcounts2(rm.spk_x,rm.spk_y,'XBinEdges',x_edge,'YBinEdges',y_edge)));
            rm.unsmoothed = spk_in_bin./occp_time;
            rm.unsmoothed(isnan(rm.unsmoothed)) = 0; 
            
            % coherence based on Fisher-z transform
            autocorr1 = lib.Rate_Matrix.cross_corr(rm.unsmoothed, rm.unsmoothed, 'stride', [bin_num - 1,bin_num - 1]);
            r = (sum(autocorr1, 'all') - 1)/8;
            rm.coherence = 0.5*log((1+r)/(1-r));
            
            % gaussian kernel
            sigma = min([delta_x, delta_y])/3;
            Gauss_K = @(v1, v2) exp(- (v1.^2 + v2.^2)/2*sigma^2);
            
            t = min(rm.time) : delta_t : max(rm.time);
            p_x = @(t) rm.itp_x(t);
            p_y = @(t) rm.itp_y(t);
                
            % spatial-smoothed rate map
            s_x = min(rm.x) + delta_x/2 : delta_x : max(rm.x) - delta_x/2;
            s_y = min(rm.y) + delta_y/2 : delta_y : max(rm.y) - delta_y/2;
            [S_x, S_y] = meshgrid(s_x, s_y);
            rm.rate_map = arrayfun(@(x, y) sum(Gauss_K(rm.spk_x - x, rm.spk_y - y))/trapz(Gauss_K(p_x(t) - x, p_y(t) - y)), S_x, S_y);
            
            % spatial information of the cell
            occp_prob = occp_time/(max(rm.time) - min(rm.time));
            rm.sp_info = sum(occp_prob.*rm.rate_map.*log2(rm.rate_map/mean(rm.rate_map,'all')),'all');
            
            % sparsity of the cell
            rm.sparsity = mean(rm.rate_map,'all')^2/sum(occp_prob.*rm.rate_map.^2, 'all');
        end
        
        % autocorrelogram of rate map
        function rm = cal_auto_corr(rm)
            rm.auto_corr_m = lib.Rate_Matrix.cross_corr(rm.rate_map, rm.rate_map, 'shape', 'same');
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
            
            % x,y centering and denoising
            x = smoothdata(x - min(x), 'gaussian', 15);
            y = smoothdata(y - min(y), 'gaussian', 15);
            
            % drop invalid time of spikes
            vld = logical(spk_time <= max(time)) & logical(spk_time >= min(time));
            spk_time = spk_time(vld);
            
            % time centering
            spk_time = spk_time - min(time);
            time = time - min(time);
        end
        
        % cross correlation matrix of 2 rate matrix
        function corr_m = cross_corr(rm1, rm2, varargin)
            defaultShape = 'full';
            defaultStride = [1,1];
            
            p = inputParser;
            addRequired(p, 'rm1');
            addRequired(p, 'rm2');
            addParameter(p, 'shape', defaultShape);
            addParameter(p, 'stride', defaultStride);
            parse(p, rm1, rm2, varargin{:});
            
            flp_rm2 = flip(flip(p.Results.rm2,2));
            corr_m = conv2(p.Results.rm1, flp_rm2, p.Results.shape)/sqrt(sum(p.Results.rm1.^2,'all')*sum(p.Results.rm2.^2,'all'));
            corr_m = corr_m(1:p.Results.stride(1):end,1:p.Results.stride(2):end);
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
        function plot_rate_map(rm, varargin)
            defautlPltUnvisited = false;
            defaultPltUnsmoothed = false;
            
            p = inputParser;
            addRequired(p, 'rm')
            addParameter(p, 'plt_unvisited', defautlPltUnvisited);
            addParameter(p, 'plt_unsmoothed', defaultPltUnsmoothed);
            parse(p, rm, varargin{:});
            
            if p.Results.plt_unsmoothed
                rate_map = p.Results.rm.unsmoothed;
            else
                rate_map = p.Results.rm.rate_map;
            end
            
            if p.Results.plt_unvisited
                rate_map(p.Results.rm.visits <= 0) = NaN;
            end
            
            figure;
            hold on;
            title(['place field of cell ',num2str(p.Results.rm.cell_id),'; SI: ', num2str(p.Results.rm.sp_info,'%.4f'), ', SP: ', num2str(p.Results.rm.sparsity, '%.4f'), ', CH: ', num2str(p.Results.rm.coherence, '%.4f')]);
            heat_map = imagesc(rate_map);
            set(gca,'YDir','normal');
            colormap('turbo');
            
            % set NaN value to be white
            set(heat_map,'AlphaData',~isnan(heat_map.CData));
            colorbar;
        end
        
        % plot the autocorrelogram of rm
        function plot_auto_corr(rm)
            figure;
            hold on;
            title(['autocorrelogram of cell ',num2str(rm.cell_id)]);
            hm=imagesc(rm.auto_corr_m);
            set(gca,'YDir','normal');
            colormap('turbo');
            set(hm,'AlphaData',~isnan(hm.CData));
            colorbar;
        end
        
    end
    
end