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
        is_pf;
        auto_corr_m;
        periodicity;
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
            rm.bin_num = bin_num;
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
            rm.unsmoothed(isinf(rm.unsmoothed)) = 0;
            
            % coherence based on Fisher-z transform
            autocorr1 = lib.Rate_Matrix.norm_cross_corr(rm.unsmoothed, rm.unsmoothed, 'window_size', [3, 3]);
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
            
            % determined place field
            rm.is_pf = lib.Rate_Matrix.bw_mask(rm.rate_map, 0.2*max(rm.rate_map,[], 'all'), 100);
            
            % spatial information of the cell
            occp_prob = occp_time/(max(rm.time) - min(rm.time));
            rm.sp_info = sum(occp_prob.*rm.rate_map.*log2(rm.rate_map/mean(rm.rate_map,'all')),'all');
            
            % sparsity of the cell
            rm.sparsity = mean(rm.rate_map,'all')^2/sum(occp_prob.*rm.rate_map.^2, 'all');
        end
        
        % autocorrelogram of rate map
        function rm = cal_auto_corr(rm)
            [m, n] = size(rm.rate_map);
            
            % set center window [m,n] to be [max odd number <= m, max odd number <= n]
            rm.auto_corr_m = lib.Rate_Matrix.norm_cross_corr(rm.rate_map, rm.rate_map, 'window_size', [m - mod(m+1,2), n - mod(n+1,2)]);
        end
        
        % spatial periodicity of autocorrelogram
        function rm = cal_periodicity(rm)
            % map autocorrelogram to actual coordinates
            [m, n] = size(rm.rate_map);
            delta_x = (max(rm.x) - min(rm.x)) / (m - mod(m+1,2));
            delta_y = (max(rm.y) - min(rm.y)) / (n - mod(n+1,2));
            prog_x = (min(rm.x) - max(rm.x) + delta_x)/2 : delta_x : (max(rm.x) - min(rm.x) - delta_x)/2;
            prog_y = (min(rm.y) - max(rm.y) + delta_y)/2 : delta_y : (max(rm.y) - min(rm.y) - delta_y)/2;
            
            % gridded interpolant function of autocorrelogram
            Z = rm.auto_corr_m;
            [X, Y] = ndgrid(prog_x, prog_y);
            F = griddedInterpolant(X, Y, Z);
            
            % divide the autocorrelogram into 360 degree
            % seperate each direction under fixed frequency
            sample_freq = 100;
            delta_rh = min([max(prog_x), max(prog_y)])/sample_freq;
            th = 0 : pi/180 : 2*pi;
            rh = 0 : delta_rh : min([max(prog_x), max(prog_y)]);
            
            % change polar grids into Cartesian grids and sampling
            [Th, Rh] = ndgrid(th, rh);
            [S_x, S_y] = pol2cart(Th, Rh);
            S_z = F(S_x, S_y);
            
            % normalized spatial periodicity
            xcov_Z = S_z * S_z';
            for i = 1:361
                xcov_Z(i,:) = circshift(xcov_Z(i,:), 1-i);
            end
            rm.periodicity = (sum(xcov_Z) - mean(S_z,'all')^2*length(th)*length(rh))/(sum(S_z.^2, 'all') - mean(S_z,'all')^2*length(th)*length(rh));
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

        % masking a binarized matrix
        function masked_m = bw_mask(matrix, threshold, least_num)
            % masking value lower than given threshold 
            masked_m = matrix >= threshold;
            
            % masking connected component less than given number
            CC = bwconncomp(masked_m, 4);
            invld_bin = cell2mat((CC.PixelIdxList(cellfun('length', CC.PixelIdxList) <= least_num)).');
            
            masked_m(invld_bin) = 0;
        end
        
        % normalized cross correlation matrix of 2 rate matrix
        function xcorr_m = norm_cross_corr(rm1, rm2, varargin)
            defaultWindowSize = size(rm1) + size(rm2) - 1;
            
            p = inputParser;
            addRequired(p, 'rm1');
            addRequired(p, 'rm2');
            addParameter(p, 'window_size', defaultWindowSize);
            parse(p, rm1, rm2, varargin{:});
            
            xcorr_m = normxcorr2(p.Results.rm1, p.Results.rm2);
            
            % choose the m-rows n-cols submatrix in the center of xcorr_m
            [rows, cols] = size(xcorr_m);
            center_rows = floor(rows / 2) - floor(p.Results.window_size(1) / 2) + (1:p.Results.window_size(1));
            center_cols = floor(cols / 2) - floor(p.Results.window_size(2) / 2) + (1:p.Results.window_size(2));
            xcorr_m = xcorr_m(center_rows, center_cols);
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
            defaultPltIsPF = false;
            
            p = inputParser;
            addRequired(p, 'rm');
            addParameter(p, 'plt_unvisited', defautlPltUnvisited);
            addParameter(p, 'plt_unsmoothed', defaultPltUnsmoothed);
            addParameter(p, 'plt_is_pf', defaultPltIsPF);
            parse(p, rm, varargin{:});
            
            if p.Results.plt_unsmoothed
                rate_map = p.Results.rm.unsmoothed;
            else
                rate_map = p.Results.rm.rate_map;
            end
            
            if p.Results.plt_is_pf
               rate_map = rate_map.*p.Results.rm.is_pf;
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
        
        % plot the periodicity of rm
        function plot_periodicity(rm)
            figure
            hold on;
            title(['periodicity of cell ',num2str(rm.cell_id)]);
            xlabel('degree');
            plot(0:360, rm.periodicity);
        end
        
    end
    
end