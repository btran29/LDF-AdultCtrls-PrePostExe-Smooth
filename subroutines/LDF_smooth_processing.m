%% LDF PORH processing (smoothed tracings)
% Smoothes laser doppler data, then finds max perfusion and fractions of
% max perfusion. Findings are output as an array if batch_processing
% condition is logical fasle. Otherwise, data is saved and exported as
% figures.
% 
% The script presumes that the data follows the model tracing for a typical
% post-occlusive reactive hyperemia test, as shown here:
% http://www.perimed-instruments.com/post-occlusive-reactive-hyperemia
% and here:
% http://www.ncbi.nlm.nih.gov/pubmed/15709655
%

%% Input
% Declare static variables necessary for smoothed LDF processing, and limit
% data to 10 minutes post-occlusion. To use this script for single studies,
% first un-comment the 'startOcc' and 'endOcc' occlusion timing variables
% and ensure they line up with the current study of interest. Second,
% change 'batch_processing' to 0.

% Batch processing mode
batch_processing = 1; 
% 1 = disable output in MATLAB; enable export; use supplied values for occlusion timings
% 0 = enable output in MATLAB; disable export; use values below for occlusion timings

% % Required when processing single studies
if batch_processing == 0
startOcc = 375.8480;	
endOcc = 436.0410;
end

% Misc
endData     = 9999; %last timepoint used to locate post-occ plateau range
rSquared    = 0.25; %for finding post-occlusion plateau (from 0-1)
delayedTime = 30; %time(sec) at which max peak is marked as excessively delayed 
maxLimit    = 30; %time(s) to limit finding max after stated end occlusion

% Smooth raw data
smPerf = smooth(perfusion,(61));

% Max 10 minutes post-occlusion
indTenminTime = find(time<endOcc+600);
tenminTime    = time(indTenminTime);
smPerf        = smPerf(indTenminTime);

%% Pre-occlusion data processing

% Select 2 minutes pre-occlusion
[idxTPreOcc] = find(tenminTime > (startOcc - 120) & tenminTime < startOcc);
tPreOcc      = tenminTime(idxTPreOcc);
pPreOcc      = smPerf(idxTPreOcc);

% Mean of selected data
plat1 = mean(pPreOcc);

%% Post-occlusion data processing
% Skip remaining processing steps if the perfusion peak is > 25 seconds
% after the occlusion as ended. Otherwise, continue to find fractions of
% max values.

% Select post-occlusion data
[idxTPostocc] = find(tenminTime > endOcc & tenminTime < endData);
tPostOcc      = tenminTime(idxTPostocc);
pPostOcc      = smPerf(idxTPostocc);

% Limit finding max perfusion within X seconds of end-occlusion
[~,idxMaxtimeRestr] = min(abs(tPostOcc - (endOcc + maxLimit))); % changed from 120

% Find max perfusion
[pMax,idxMax] = max(pPostOcc(1:idxMaxtimeRestr));
tMax          = tPostOcc(idxMax);

% Is the peak 'delayed?'
% Mark study if pMax occurs (delayedTime) seconds after endOcc
if abs(tPostOcc(idxMax) - endOcc) >= delayedTime;
    delayed_peak = 1;
else
    delayed_peak = 0;
end

%% Locate post-occlusion plateau with linear interpolation
% Find the timepoint at which the post-occlusive decay stops. The current
% method involves loosly fitting a third-degree polynomial to the decay
% area followed by locating the first point at which the second derivative
% of that polynomial passes 0 (i.e. when the decay stops or the inflection
% point of the first concavity in the decay function).
%

for i = 3:8 % loosely fit poly via SLMengine to post-occ to r_squared
    slm_plat2 = slmengine(tPostOcc(idxMax:end),pPostOcc(idxMax:end),'degree',3,...
        'knots',i,'interior','free','plot','off');
    if slm_plat2.stats.R2 >= rSquared %refer to slmengine doc
        break;
    end
end

% Find 2nd deriv of fit curve (if hyperemia isn't monotonic)
slmPlat2d2 = slmeval(tPostOcc, slm_plat2,2);

% Locate plateau 2
[~,idxPlat2Start]  = min(abs(slmPlat2d2));
tPlat2Start        = tPostOcc(idxPlat2Start);
idxPlat2           = tPostOcc > tPostOcc(idxPlat2Start) & ...
                      tPostOcc < (tPostOcc(idxPlat2Start) + 120);
                  
% Mean of post-occ plateau perf
plat2 = mean(pPostOcc(idxPlat2));


% Find minimum of 1st deriv of fit curve (if hyperemia is monotonic)
if isnan(plat2) == 1
    % Find 1st derivative
    slmPlat2d1 = slmeval(tPostOcc, slm_plat2,1);
    
    % Locate plateau 2
    [~,idxPlat2Start] = min(abs(slmPlat2d1));
    
    % Set plateau 2 to 2 min in length
    idxPlat2   = tPostOcc > tPostOcc(idxPlat2Start) & ...
                      tPostOcc < (tPostOcc(idxPlat2Start) + 120);
                  
    % Mean of post-occ plateau perf
    plat2 = mean(pPostOcc(idxPlat2));
    
    % Obtain plat2 time
    tPlat2Start = tPostOcc(idxPlat2Start);
end

% Find area under the curve from occlusion to plat 2 

% Adjust data via proxy variable to simplify function call
pBloodFlow = pPostOcc;

% Remove all values below plat2 for AUC
pBloodFlow(pBloodFlow <= plat2) = plat2;

% Normalize data to plat2
pBloodFlow = pBloodFlow - plat2;

% Integral functions only accept function handles at outputs, thus use
% trapz function to generate trapezoids for the AUC calc
bloodflow = trapz(tPostOcc(idxMax:idxPlat2Start),pBloodFlow(idxMax:idxPlat2Start));

% Skip further analysis if plateau 2 is greater than max
if plat2 > pMax
    delayed_peak = 1; %#ok<NASGU>
    continue  %#ok<CONT>
end

%% Locating fractions of max flow with linear interpolation
% Locate points directly before and after the 75, 50, and 25% of max in the
% perfusion dataset. If fractions of max are close together, try looking at
% the 2nd plateau value: it may be close to pMax.

maxFrac   = [75 50 25];
tMaxFracs = cell(size(maxFrac));
pMaxFracs = cell(size(maxFrac));

% restrict finding fractions of max after max time point
maxFracPerfPostOcc = pPostOcc(idxMax:end);
maxFracTimePostOcc = tPostOcc(idxMax:end);

for iMaxFrac = 1:numel(maxFrac)
    % the current target fraction to solve for
    target = pMax-(pMax-plat2)*(1-(maxFrac(iMaxFrac)/100));

    % simplify function by reducing target level
    maxFracPerfPostOccTemp = maxFracPerfPostOcc - target;

    % find crossings
    indInterp = find(maxFracPerfPostOccTemp(1:end-1).*maxFracPerfPostOccTemp(2:end)<0);

    % calculate interpolated X values
    Xinterp = (maxFracTimePostOcc(indInterp) + (maxFracTimePostOcc(indInterp+1) - maxFracTimePostOcc(indInterp))...
        .* maxFracPerfPostOccTemp(indInterp) ./ (maxFracPerfPostOccTemp(indInterp) - maxFracPerfPostOccTemp(indInterp+1)));

    % find closest point in dataset that is closest to first Xinterp
    [~,indMaxFrac] = min(abs(Xinterp(1)-maxFracTimePostOcc));

    % insert new data
    pMaxFracs{iMaxFrac} = maxFracPerfPostOcc(indMaxFrac);
    tMaxFracs{iMaxFrac} = maxFracTimePostOcc(indMaxFrac);
end

% Time values associated with fractions of max flow
tTQuartermax = tMaxFracs{1};
tHalfmax     = tMaxFracs{2};
tQuartermax  = tMaxFracs{3};

% Perfusion values associated with fractions of max  flow
tquartermax = pMaxFracs{1};
halfmax     = pMaxFracs{2};
quartermax  = pMaxFracs{3};

%% Output

% Located variables

% Temperature
%     temp

% Perfusion
%     plat1
%     plat2
%     pMax
%     tquartermax
%     halfmax
%     quartermax
%     bloodflow

% Time
%     tMax
%     tTquartermax
%     tHalfmax
%     tQuartermax
%     tPlat2Start

% Calculated variables
percentDiffPlats      = ((plat2-plat1)/(plat1))*100;
plat2OvrPlat1         = plat2/plat1;
pMaxOvrPlat2          = pMax/plat2;
pMaxOvrPlat1          = pMax/plat1;
endOccTo_tMax         = tMax - endOcc;
endOccTo_tTQuartermax = tTQuartermax - endOcc;
endOccTo_halfmax      = tHalfmax - endOcc;
endOccTo_tQuarterMax  = tQuartermax - endOcc;
tMaxTo_tTQuartermax   = tTQuartermax - tMax;
tMaxTo_halfmax        = tHalfmax - tMax;
tMaxTo_tQuarterMax    = tQuartermax - tMax;
perfVel               = abs(pMax-plat1)/(tMax - endOcc);
stdPostOcc            = std(pPostOcc(idxPlat2));
avgTemp               = mean(temp);

% Figures
% Generate troubleshooting plot and SLM plot, used to determine the
% post-occlusive plateau. In single-study mode, pop-up associated figures. 
% In batch-processing mode, disable pop-up and automatically save figures
% into new directories.


close all

% Temporary time indicies for quicker plotting
rf2g(1:numel(tPostOcc))   = plat2;
rf1g(1:numel(idxTPreOcc)) = plat1;
twominuteindex_set        = find(tenminTime>endOcc & tenminTime<endOcc+30);
avgTempg(1:numel(tenminTime)) = avgTemp;

% Current directory
currDir = pwd; % make folders + set current Dir if they don't exist
if batch_processing == 1 && exist([pwd '\Plat2 SLM Plots'],'dir') == 0
    mkdir('Plat2 SLM Plots')
    mkdir('Troubleshooting Plots')
end

% SLM Plot
plotslm(slm_plat2);
title(sprintf('PostOcc fit for %s',study));

% Export SLM plot
if batch_processing == 1;
    set(gcf,'Visible','off', 'Color', 'w');
    cd([currDir '\Plat2 SLM Plots'])
    export_fig(sprintf('PostOcc fit for %s',study),'-png','-m2');
    cd(currDir)
else
    set(gcf,'Visible','on', 'Color', 'w');
end

% Troubleshooting Plot
figure('name', sprintf('Troubleshoot overview for %s',study));

% Section Markers subplot
subplot(2,1,1);
hold on;

% Data
plot(tenminTime,smPerf,...
     tenminTime,pressure(indTenminTime),'m-',...
     tenminTime,temp(indTenminTime),...
     tPostOcc,pPostOcc,'r-',...
     tenminTime(idxTPreOcc),rf1g(1:length(tenminTime(idxTPreOcc))),'k--',...
     tPostOcc(1:end),rf2g(1:length(tPostOcc)),'k--');
 
% Labels
xlabel('Time (seconds)');
ylabel('Perfusion Units (PU)');
legend('Perf','Press','Temp','PostOccPerf');
title(sprintf('Data selections for %s',study));

% Fractions of max zoomed-in view subplot
subplot(2,1,2);
hold on;

% Fractions of max
scatter(tMax,pMax,75,'k','markerfacecolor','k');
scatter(tTQuartermax,tquartermax,75,'dk','markerfacecolor','k');
scatter(tHalfmax,halfmax,75,'vk','markerfacecolor','k');
scatter(tQuartermax,quartermax,75,'sk','markerfacecolor','k');

% Post-occ data
plot(tenminTime(twominuteindex_set),smPerf(twominuteindex_set),'b-',...
    tenminTime(twominuteindex_set),...
    rf2g(1:length(tenminTime(twominuteindex_set))),...
    'k--','MarkerSize',12);

% Labels
xlabel('Time (seconds)');
ylabel('Perfusion Units (PU)');
legend('Max','75%max','50%max','25%max','Perf','PostOccPlat');

% Re-shift bloodflow AUC
pBloodFlowGraph = pBloodFlow + plat2;

% Initialize verticies
xverts = [transpose(tPostOcc(idxMax:idxPlat2Start-1));transpose(tPostOcc(idxMax:idxPlat2Start-1));...
          transpose(tPostOcc(1+idxMax:idxPlat2Start));transpose(tPostOcc(1+idxMax:idxPlat2Start))];
yverts = [repmat(plat2,1,numel(idxMax:idxPlat2Start-1));transpose(pBloodFlowGraph(idxMax:idxPlat2Start-1));...
          transpose(pBloodFlowGraph(1+idxMax:idxPlat2Start));repmat(plat2,1,numel(idxMax:idxPlat2Start-1))];

% Draw trapezoids
p = patch(xverts,yverts,'b','EdgeAlpha',.00,'FaceAlpha',0.1);
axis ([tenminTime(twominuteindex_set(1)) tenminTime(twominuteindex_set(end)) 0 inf])

% Export
if batch_processing == 1;
    set(gcf,'Visible','off', 'Color', 'w');
    cd([currDir '\Troubleshooting Plots'])
    export_fig(sprintf('Troubleshoot for %s',study),'-png','-m2');
    cd(currDir)
else
    set(gcf,'Visible','on', 'Color', 'w');
end