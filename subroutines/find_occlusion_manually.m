%%Find Occlusion Mannually
% Creates  a figure to help find occlusion timings manually.
% Assumes conversion from csv to mat; only works on imported 
% variables in mat .mat files
close all
figure
set(gcf,'Visible','on', 'Color', 'w');
plot(time,perfusion,time,pressure)
title(sprintf('%s',study))
axis([-inf inf 0 inf])
xlabel('Time (seconds)');
ylabel('Perfusion Units (PU)');
fprintf('Plotting %s \n',study); 
clear all

%% Locating pre and post occlusion using pressure data
% Need to optimize for speed; inaccurate by a few milliseconds. Uses shape 
% language modeling tool to find 'breakpoints' of a piecewise linear fit 
% curve to given pressure data.
%
%for i=length(matfiles.name);
% occRSquared = .95;
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
%     slmOccD1 = slmeval(time, slm_occ,1); %first derivative of fit curve
%     [~,idxStartOcc]= max(slmOccD1);
%     [~,idxEndOcc] = min(slmOccD1);   
%     startOcc = time(idxStartOcc);
%     endOcc = time(idxEndOcc);
%end