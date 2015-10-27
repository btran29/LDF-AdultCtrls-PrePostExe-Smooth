%Figure
figure('Name',sprintf('Overview for %s\n',study)); hold on;
set(gcf,'Visible','off', 'Color', 'w');
plot(time,perfusion,time,pressure);
xlabel('Time (seconds)');
ylabel('Perfusion Units (PU)');
legend('Perfusion','Pressure');

% Make new folder to save tracings if none exists
currDir = pwd;
if batch_processing == 1 && exist([pwd '\Processed LDF Tracings'],'dir') == 0
    mkdir('Processed LDF Tracings') 
end 

% Save to new folder with batch processing
if batch_processing == 1 
    cd([currDir '\Processed LDF Tracings'])
    export_fig(sprintf('Processed LDF Tracing for %s',study),'-png','-m2');
    cd(currDir)
else % Save to current folder w/o batch processing
    export_fig(sprintf('Processed LDF Tracing %s',study),'-png','-m2');
end

