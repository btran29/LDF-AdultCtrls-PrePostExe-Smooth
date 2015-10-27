% Laser Doppler Processing for Peaks only
start_occ = 196.727;
end_occ = 257.447;
minpkdis = 14; %sets minimum peak distance for rest flow and post-occlusion peak analysis, best to get low enough to get all major peaks
time_restr = 9999; %last timepoint to find max and fractions of max
perf_limit = 1000; % in case you need to remove high PU noisy peaks
r_squared = 0.25; %for finding post-occlusion plateau (from 0-1)
end_data = 9999; %last timepoint used to locate post-occ plateau range


batch_processing = 0; %for batch processing 1 = disable output for ea study
hampel = 0; %enable hampel filtering @ 10 sec window, 5std dev
%% Locating pre and post occlusion using pressure data
% occ_r_squared = .95;
%     for i = 3:20; %repeats slm operation until the curve fits above a specified r-squared
%         fprintf('\t Locating pre/post occlusion with i = %g\n',i)
%         slm_occ = slmengine(time, pressure, 'degree', 1, 'knots', i, 'interior', 'free', 'plot', 'off');
%         if slm_occ.stats.R2>= occ_r_squared; %r-squared threshold
%             sprintf('\t found pre/post occ successfully')
%             break
%         end
%         if i==20;
%             sprintf('\t hit r_squared ceiling; consider finding start/end occ manually \n')
%             break
%         end
%     end
%     slm_occ_d1 = slmeval(time, slm_occ,1); %first derivative of fit curve
%     [~,idx_startocc]= max(slm_occ_d1);
%     [~,idx_endocc] = min(slm_occ_d1);   
%     start_occ = time(idx_startocc);
%     end_occ = time(idx_endocc);  
%% Pre-occlusion data processing
% Selection of timepoints
[pre_occ_time_idx]=find(time>(start_occ-120) & time<start_occ); %2 minutes of pre-occlusion data
pre_occ_time = time(pre_occ_time_idx);
pre_occ_perf = perfusion(pre_occ_time_idx);

% Excluding laserdata above threshold
[pre_occ_outliers] = find(pre_occ_perf>perf_limit);
pre_occ_perf(pre_occ_outliers) = NaN;
%
% Peak processing
[~,pre_occ_perf_pks_idx]=findpeaks(pre_occ_perf,'minpeakdistance',minpkdis);
rf1 = mean(pre_occ_perf(pre_occ_perf_pks_idx));
%% Post-occlusion data processing
% Selection of time points
[post_occ_time_idx]= find(time>end_occ & time<end_data); %begin/endpoint
time_set = time(post_occ_time_idx);
perf_set = perfusion(post_occ_time_idx); 

% Removing sensor movement artifacts via hampel filter
% [YY,I,Y0,LB,UB] = hampel(time_set,perf_set,5,5);figure; %10-sec win;5std dev

% Excluding laserdata above threshold
[post_occ_outliers] = find(perf_set>perf_limit);
perf_set(post_occ_outliers) = NaN;

% Selecting data from max perfusion to end_data
[~,idx_maxtime_restr] = min(abs(time_set-(time_restr))); %limits finding max to time-restr
if time_set(idx_maxtime_restr) > (end_occ+30) %limits finding max to 30 sec post end_occ
    [~,idx_maxtime_restr] = min(abs(time_set-(end_occ+30))); 
end
[~,idx_max] = max(perf_set(1:idx_maxtime_restr));
[newtime_i3]= find(time_set>=time_set(idx_max)-1); %include 1 seconds before max for circshift
time_set = time_set(newtime_i3);
perf_set = perf_set(newtime_i3);
%% Finding fractions of max flow
% Peak processing
[peaks,location]=findpeaks(perf_set,'minpeakdistance',minpkdis); % plot(time_set(location),perf_set(location),'m')
time_peaks = time(location);
perf_peaks = perf_set(location);

% Locating post-occlusion plateau
for i = 3:10; %repeats slm operation until the curve fits above a specified r-squared
    slm_plat2 = slmengine(time_peaks, perf_peaks, 'degree', 3, 'knots', i, 'interior', 'free', 'plot', 'off');
    if slm_plat2.stats.R2 >= r_squared; %r-squared threshold
        break;
    end;
end;
slm_plat2d2 = slmeval(time_peaks, slm_plat2,2); %second derivative of fit curve

% select points before and after passing through 0
after    = (circshift(slm_plat2d2<0,1) & slm_plat2d2>0) | (circshift(slm_plat2d2>0,1) & slm_plat2d2<0);
after(1) = false;
before   = circshift(after,-1);
xx = [time_peaks(before) time_peaks(after)];
yy = [slm_plat2d2(before) slm_plat2d2(after)];

% interpolate and insert new data
plat2_start = time_peaks(before) - (slm_plat2d2(before)).*diff(xx,[],2)./diff(yy,[],2); 
if isempty(plat2_start) == 1;
    [~,plat2_start_idx]= min(abs(slm_plat2d2));
    plat2_idx = find(time_peaks>time_peaks(plat2_start_idx) & time_peaks<(time_peaks(plat2_start_idx)+120));
else
    plat2_idx = find(time_peaks>plat2_start(1) & time_peaks<(plat2_start(1)+120));
end
if isempty(plat2_idx) == 1; %use last 2 minutes of previous methods fail
    plat2_idx = find(time_peaks>time_peaks(end)-120);
end
rf2 = mean(perf_peaks(plat2_idx));

% Max-related variables
[~,idx_maxtime_restr] = min(abs(time_peaks-(end_occ+60))); %limits finding max to 60 sec after end_occ
[perf_max,idx_max] = max(perf_set(1:idx_maxtime_restr));

% delayed peak identification
if abs(time_set(idx_max)-(end_occ)) >= 15;
    delayed_peak = 1;
else
    delayed_peak = 0;
end

% y == 75%, 50%, 25% limits
lims = [75 50 25];
X = cell(size(lims));
Y = cell(size(lims));
for ind = 1:numel(lims)

    % the current level limit to solve for
    lim = perf_max-(perf_max-rf2)*(1-(lims(ind)/100));

    % select points before and after passing through the current limit
    after    = (circshift(perf_peaks<lim,1) & perf_peaks>lim) | (circshift(perf_peaks>lim,1) & perf_peaks<lim);
    after(1) = false;
    before   = circshift(after,-1);
    xx = [time_peaks(before) time_peaks(after)];
    yy = [perf_peaks(before) perf_peaks(after)];

    % interpolate and insert new data
    new_X = time_peaks(before) - (perf_peaks(before)-lim).*diff(xx,[],2)./diff(yy,[],2);
    X{ind} = new_X;
    Y{ind} = lim * ones(size(new_X));

end

%Max variables
time_max_offset = time_set(idx_max)-time_set(1);
time_idx_tquartermax = find(X{1}>(time_max_offset),1,'first');
time_idx_halfmax = find(X{2}>(time_max_offset),1,'first');
time_idx_quartermax = find(X{3}>(time_max_offset),1,'first');

tquartermax = Y{1}(time_idx_tquartermax);
halfmax = Y{2}(time_idx_halfmax);
quartermax = Y{3}(time_idx_quartermax);

time_tquartermax = X{1}(time_idx_tquartermax)+time_set(1);
time_halfmax = X{2}(time_idx_halfmax)+time_set(1);
time_quartermax = X{3}(time_idx_quartermax)+time_set(1);

% %Means of first and last two minutes of data
% mean1_idx = find(pre_occ_time<=(pre_occ_time(1)+120));
% mean1 = mean(pre_occ_perf(mean1_idx));
% mean2_idx = find(time_set>=(time_set(end)-120));
% mean2 = mean(perf_set(mean2_idx));
%% Output
close all
% Temporary variables for quicker plotting
rf2g(1:length(time_set)) = rf2; 
rf1g(1:length(pre_occ_time_idx)) = rf1;
twominuteindex_set = find(time_set<(time_quartermax+10));
twominuteindex_pks = find(time_peaks<(time_quartermax+10-time_set(1)));

%Data in lieu of batch processing (may be outdated; check with current data
%table!) updated: 6-9-14
if batch_processing == 0;
        outputVars = [studyNameOutput(1) ...
                      studyNameOutput(2) ...
                      studyNameOutput(3) ...
                      start_occ ...
                      end_occ ...
                      time_restr ...
                      minpkdis ...
                      r_squared ...
                      perf_limit ...
                      end_data ...
                      rf1 ...
                      rf2 ...
                      perf_max ...
                      time_set(idx_max) ...
                      tquartermax ...
                      time_tquartermax ...                      
                      halfmax ...  
                      time_halfmax ...                      
                      quartermax ...                      
                      time_quartermax ...
                      (rf2/rf1) ...
                      (perf_max/rf2) ...
                      (perf_max/rf1) ...
                      (time_set(idx_max) - end_occ) ... % Time Variables
                      (time_tquartermax - end_occ) ...
                      (time_halfmax - end_occ) ...
                      (time_quartermax - end_occ) ...
                      (time_tquartermax - time_set(idx_max)) ...
                      (time_halfmax - time_set(idx_max)) ...
                      (time_quartermax - time_set(idx_max)) ... % Perfusion Variables
                      ];
end
%SLM Plot
plotslm(slm_plat2);
    title(sprintf('PostOcc fit for %s',study));
    if batch_processing == 1;
        set(gcf,'Visible','off', 'Color', 'w');
        export_fig(sprintf('02 PostOcc fit for %s',study),'-png','-m2');
    else 
        set(gcf,'Visible','on', 'Color', 'w');
    end
%Troubleshoot plot
figure('name', sprintf('Troubleshoot overview for %s\n',study));
    subplot(2,1,1);plot(time,perfusion,...
                         time,pressure,'m-',...
                         pre_occ_time(pre_occ_perf_pks_idx),pre_occ_perf(pre_occ_perf_pks_idx),...
                         time_set,perf_set,...
                         time(pre_occ_time_idx),rf1g,'k--',...
                         time_set,rf2g,'k--',...
                         (time_peaks(plat2_idx)+time_set(1)),perf_peaks(plat2_idx));
        hold on
        title(sprintf('Data selections for %s',study));
        if exist('idx_startocc','var')==1
             scatter(start_occ,pressure(idx_startocc));
             scatter(end_occ,pressure(idx_endocc));
        end
    subplot(2,1,2); 
        hold on;
        scatter(time_set(idx_max),perf_set(idx_max),75,'k','markerfacecolor','k');
        scatter(time_tquartermax,tquartermax,75,'dk','markerfacecolor','k');
        scatter(time_halfmax,halfmax,75,'vk','markerfacecolor','k');
        scatter(time_quartermax,quartermax,75,'sk','markerfacecolor','k');
        plot(time_set(twominuteindex_set),perf_set(twominuteindex_set),'b-',...
                time_peaks(1:twominuteindex_pks(end))+time_set(1),perf_peaks(1:twominuteindex_pks(end)),'r-',...
                time_set(twominuteindex_set),rf2g(twominuteindex_set),'k--','MarkerSize',12);
        xlabel('Time (seconds)');
        ylabel('Perfusion Units (PU)');
        legend('Max','75%max','50%max','25%max','Original Data','Interpolated Peak Function','PostOcc Plateau');
    if batch_processing == 1;
        set(gcf,'Visible','off', 'Color', 'w');
        export_fig(sprintf('01 Troubleshoot for %s',study),'-png','-m2');
    else
        set(gcf,'Visible','on', 'Color', 'w');
    end
%Analysis Plot
if batch_processing == 1;
    %do not output analysis plot
else
figure('Name',sprintf('Analysis for %s\n',study)); hold on;
    set(gcf,'Visible','on', 'Color', 'w');
    title(sprintf('Analysis for %s\n',study));
    scatter(time_set(idx_max),perf_set(idx_max),75,'k','markerfacecolor','k');
    scatter(time_tquartermax,tquartermax,75,'dk','markerfacecolor','k');
    scatter(time_halfmax,halfmax,75,'vk','markerfacecolor','k');
    scatter(time_quartermax,quartermax,75,'sk','markerfacecolor','k');
    plot(time_set,perf_set,'b-',time_peaks+time_set(1),perf_peaks,'r-',time_set,rf2g,'k--','MarkerSize',12);
    xlabel('Time (seconds)');
    ylabel('Perfusion Units (PU)');
    legend('Max','75%max','50%max','25%max','Original Data','Interpolated Peak Function','PostOcc Plateau');
end