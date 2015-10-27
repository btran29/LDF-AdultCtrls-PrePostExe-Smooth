%Plot all studies data without correction
matfiles = dir('*.mat');

% Make new folder to save tracings
currDir = pwd;
if exist([pwd '\Raw LDF Tracings'],'dir') == 0
    mkdir('Raw LDF Tracings')
end

for i = 1:length(matfiles)
    % Load file
    load(matfiles(i).name);
    % Figure
    set(gcf,'Visible','off', 'Color', 'w');
    plot(time,perfusion,time,pressure,time,temp)
    title(sprintf('%s',study))
    axis([-inf inf 0 inf])
    xlabel('Time (seconds)');
    ylabel('Perfusion Units (PU)');
    % Export to png
    cd([currDir '\Raw LDF Tracings'])
    export_fig(sprintf('%s',study),'-png','-m2');
    cd(currDir)
end