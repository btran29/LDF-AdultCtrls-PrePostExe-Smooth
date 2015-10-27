%% Batch Processing Script for Laser Doppler Flowmetry PORH Data
% This script runs processing scripts for all .mat files in the current 
% directory and exports named variables into output.xlsx. Data is sorted by
% total output and individual subject data, based on initials.
%
% Steps:
%
% # Complete prior formatting of LDF data into .mat files
% # Generate input.xlsx, and input occlusion timings manually for ea. study
% # Import input.xlsx (occ timings with inputs ordered by matfiles)
% # Navigate to mat files folder
% # Run and manually adjust studies as needed with output
%
% Changes to variable order must be noted and replicated in 3 sections:
%
% # Script loop
% # Variables and associated output array columns
% # Headers for excel output
%
% *WARNING*: Following completion, this script clears all workspace
% variables unrelated to the output.
%

%% Variables of Interest
% First, note the variables of interest. The variable names are taken from
% the processing scripts. If you intend to change variable order or add 
% additional variables, double check if the changes work with single
% studies. The spelling of the variables here will also be used as headers
% in the output. The variables of interest are formatted into a cell of
% strings as they are used as headers and variable descriptions as well. 
% The eval function will be used to locate the actual numbers associated 
% with the variables.

outputEvalVars = {'initial' , 'initial'; 
              'prePost' , 'Pre or post';
              'date' , 'date';
              'startOcc' , 'start of occlusion (sec)';
              'endOcc' , 'end of occlusion (sec)';
              'rSquared' , 'fit for determining post-occ plateau';
              'endData' , 'cutoff for data to be shown';
              'plat1' , 'absolute perfusion variables (plat 1)';
              'plat2' , 'mean of post-occ plateau';
              'percentDiffPlats' , 'percentage diff btwn plat';
              'pMax' , 'max perfusion';
              'tquartermax' , '75% of max perfusion';
              'halfmax' , '50% of max perfusion';
              'quartermax' , '25% of max perfusion';
              'bloodflow' , 'AUC from max to beginning of steady-state';
              'delayed_peak' , '% is the hyperemic peak >15s postocc?';
              'tMax' , 'absolute time variables';
              'tTQuartermax' , 'time of 75% of max perfusion';                     
              'tHalfmax' , 'time of 50% of max perfusion';                     
              'tQuartermax' , 'time of 25% of max perfusion';
              'tPlat2Start' , 'time when steady state is reached';
              'plat2OvrPlat1' , 'relative perfusion variables';
              'pMaxOvrPlat2' , 'ratio max perf/post-occ plat mean';
              'pMaxOvrPlat1' , 'ratio max perf/pre-occ plat mean';
              'endOccTo_tMax' , 'relative time variables';
              'endOccTo_tTQuartermax' , 'time end-occ to 75% max';
              'endOccTo_halfmax' , 'time end-occ to 50% max';
              'endOccTo_tQuarterMax' , 'time end-occ to 25% max';
              'tMaxTo_tTQuartermax' , 'time max to 75% max';
              'tMaxTo_halfmax' , 'time max to 50% max';
              'tMaxTo_tQuarterMax' , 'time max to 25% max';
              'perfVel' , 'perf velocity';
              'stdPostOcc' , 'st.dev post-occ plat from post-occ plat mean';
              'avgTemp', 'average temperature over whole test'
              };
numvars = size(outputEvalVars,1);

%% Begin Error Logging
% Saves the subsequent matlab command window text to a text file. This is 
% useful for going over studies that didn't process, or looked strange 
% following the analysis. 

% Make folder + set error log directory if it doesn't exist
currDir = pwd; 
if exist([pwd '\Error Logs'],'dir') == 0
    mkdir('Error Logs')
end

% Use current time for error log
currTime = datestr(now);
currTime = strrep(currTime,':','.');
diary ([currDir,'\Error Logs','\Error Log ',currTime,' (DateTime)','.txt'])

%% Load Input.xlsx
% Loads input timings from 'generate_input_for_batch' into a strcuture 
% variable. This assumes that the input file is located in the working 
% directory (along with the .mat files), and that the first sheet of the 
% input.xlsx file is named 'input'.

inputOccTimings = importdata('input.xlsx');
inStartOcc = inputOccTimings.data.input(:,1);
inEndOcc = inputOccTimings.data.input(:,2);


%% Batch loop
% Pre-allocate the output array based on .mat files in the 
% working directory. Then, ask that MATLAB run the specified scripts 
% on every individual .mat file and save the output into a combined array.

matfiles = dir('*.mat'); %removed minpkdis/perf_limit for smoothed data
outputArray = cell((length(matfiles)+1),numvars);

% Run specified scripts in a loop:
for iM = 1:(numel(matfiles))
try              
    load(matfiles(iM).name);
    %Load input variables for particular study
        startOcc = inStartOcc(iM);
        endOcc = inEndOcc(iM);
        %minpkdis = in_minpkdis(iM);

    %Scripts to be run
        LDF_smooth_processing;
        overview_plot;
        %LDF_peaks_processing;
         
    % Use Eval function so variable names can serve as headers as well
%     outputVars = cell(1,size(outputEvalVars,2));
%     for f = 1:length(outputEvalVars)
%         outputVars{f} = eval(outputEvalVars{f});
%     end
    outputVars = cell(1,size(outputEvalVars,1));
    for f = 1:size(outputEvalVars,1)
        outputVars{f} = eval(outputEvalVars{f,1});
    end
    outputArray(iM+1,1:size(outputEvalVars,1)) = outputVars;

    %Progress
        fprintf('(%g/%g) Processed %s \n',iM,numel(matfiles),study);  
        close all;
        clearvars -except matfiles...
                          in_time_restr...
                          inStartOcc...
                          inEndOcc...
                          in_minpkdis...
                          in_perf_limit...
                          err...
                          outputArray...
                          outputEvalVars...
                          numvars...
                          inputOccTimings
catch BatchError
    fprintf('Error on (%g/%g) %s \n',iM,numel(matfiles),study);
  continue      
end
end

%% Make a spreadsheet and save output to mat
% Write output onto excel sheet, with headers. The output variables double
% as headers, so the columns between headers and variables should always
% match.

fprintf('Starting outputArray processing...\n');
filename = 'output.xlsx';

warning('off','MATLAB:xlswrite:AddSheet');

% Headers
headers = transpose(outputEvalVars(:,1));

% Sheet Locations (individual subejct data is after grouped output)
generation_time = 1;
original_output = 2;
pruned_output = 3;

% Export current time
currTime = strsplit(datestr(now));
currTime = ['Generated', currTime(1), currTime(2)];

varsUsed = {'Variables used:',''};
varsUsed = [varsUsed; outputEvalVars];
xlswrite(filename,currTime,generation_time,'A1');
xlswrite(filename,varsUsed,generation_time,'A4');

% Move output data to another variable
grpData = outputArray;
xlswrite(filename,outputArray,original_output,'A2');
xlswrite(filename,headers,original_output,'A1');
 
%% Match variables to columns for study pruning & calculations
% Tell Matlab which columns in the output array correspond to a particular 
% variable. They MUST match in order for the pre & post study pairing, 
% post/pre calculations, and filtered mean analysis to function. The
% following code allows you change the order of the variables. The spelling
% of the variables must match the variables of interest in
% 'outputEvalVars.' An alternative to this code would be to statically
% note a variable with its corresponding column in the output array (e.g.
% colDates = 3).

calcColsHeaders = {'date', ...
                   'plat1', ... 
                   'pMax', ...
                   'pMaxOvrPlat1', ...
                   'endOccTo_tMax',...
                   'perfVel',...
                   'avgTemp',...
                   'prePost'
                   };

% Given variables of interest above, find location
calcCols = zeros(1,numel(calcColsHeaders));
for iCalcCols = 1:length(calcColsHeaders)
    str = strncmp(outputEvalVars(:,1),calcColsHeaders(iCalcCols),...
            numel(calcColsHeaders{iCalcCols}));
    calcCols(iCalcCols) = find(str,1,'first');
end

colDates        = calcCols(1); 
colPlat1        = calcCols(2);
colPMax         = calcCols(3);
colPMaxOvrPlat1 = calcCols(4);
colEndOccToMax  = calcCols(5);
colPVel         = calcCols(6);
colAvgTemp      = calcCols(7);
colPrePost      = calcCols(8);

%% Break data into groups based on unique initials & remove extraneous studies
% For further calculations, studies that do not have a corresponding
% pre or post match on the same day (and for the same subject) are removed.
% Subjects with only one study (in total) are also removed. The remaining
% studies are saved into a new 'total' data table on a new sheet and 
% grouped by subject initals on separate sheets.
%
% *Requires study date column to correspond with study date variable
%

fprintf('Pairing pre-post studies...\n');

% Remove empty rows from raw data
grpData(any(cellfun(@isempty,grpData)'),:) = [];

% Get unique initials
[uniqueInit,initStart,~] = unique(grpData(:,1),'stable');

% Pre-allocate for loop
outGrp = cell(numel(uniqueInit),1);

% Loop for N-1 subjects
for iSubj = 1:length(initStart)-1
    startGrp = initStart(iSubj);
    endGrp   = initStart(iSubj+1);
    currGrp  = grpData(startGrp:endGrp-1,:);
    outGrp(iSubj,1) = {currGrp};
end
    startGrp = initStart(end);
    endGrp   = size(grpData,1);
    currGrp  = grpData(startGrp:endGrp,:);
    outGrp(length(initStart),1) = {currGrp};
    
% Mark studies with missing data or EndOccToMax >15 sec
removeCell{1,1} = true; % mark cells with logical 'true'
for iSubj = 1:size(outGrp,1);
    for iRow = 1:size(outGrp{iSubj,1},1) 
        currState = cell2mat(outGrp{iSubj,1}(iRow,colPrePost));
        
        % Remove studies with unknown pre or post-exe status
        if strncmp(currState,'Unk',3) == 1  
            outGrp{iSubj,1}(iRow,:) = removeCell;
        end
        
        % Remove studies with EndOccToMax > 15 from pruned data
        currEndOccToMax = cell2mat(outGrp{iSubj,1}(iRow,colEndOccToMax));
        if currEndOccToMax > 15
            outGrp{iSubj,1}(iRow,:) = removeCell;
        end
    end
    % Remove marked rows    
    outGrp{iSubj,1}(any(cellfun(@islogical,outGrp{iSubj,1})'),:) = [];
end

% After removing studies with missing data, find studies without a
% corresponding pre or post-exercise state
clearvars duplDates
duplDates = cell(size(outGrp,1),1);

for iSubj = 1:size(outGrp,1) %Subject Loop
    dates = outGrp{iSubj,1}(:,colDates);
    % Find number of occurances for each date
    [~,~,dateOccur] = unique(dates); 

    % Initialize array to mark duplicates
    z = size(dateOccur,1); 

    % Mark Pairs
    for iDates = 1:length(dateOccur)
        z(iDates) = sum(dateOccur == dateOccur(iDates));
    end
    
    for iDates = 1:length(z);
        if z(iDates) ~=2 % mark studies w/ 1 or >2 occurrences on the same date
            duplDates{iSubj,1}(iDates,1) = 1;
        else
            duplDates{iSubj,1}(iDates,1) = 0;
        end
    end % end marking loop
end % end subject loop

% Mark rows in output array
removeCell{1,1} = true; % mark cells with logical 'true'
for iSubj = 1:size(outGrp,1);
    for iDuplDates = 1:size(duplDates{iSubj,1},1)
        if duplDates{iSubj,1}(iDuplDates,1) == 1
            outGrp{iSubj,1}(iDuplDates,:) = removeCell;
        end
    end
end

% Remove marked rows
for iSubj = 1:size(outGrp,1)
    outGrp{iSubj,1}(any(cellfun(@islogical,outGrp{iSubj,1})'),:) = [];
end

% Remove subjects w/ 1 study (e.g. only one post-exe in total)
for iSubj = 1:size(outGrp,1)
    try
        if size(outGrp{iSubj},1) == 0;
            outGrp(iSubj) = [];
        end
        % Note removal of subjects with single studies
        catch subjRemoval
    end
end

% Compile pruned output into totals table
% Pre-allocate
matfiles = dir('*.mat');
allSubj = cell(numel(matfiles),numvars);
totSubjRows = zeros(1,1);

% Recursive Subject Loop
for iSubj = 1:size(outGrp,1)
    currSubjRows = size(outGrp{iSubj,1},1);
    allSubj(totSubjRows+1:totSubjRows+currSubjRows,:) = outGrp{iSubj,1};
    totSubjRows = totSubjRows + currSubjRows;
end

% Remove empty rows
outGrp(all(cellfun(@isempty,outGrp),2), : ) = [];
allSubj(all(cellfun(@isempty,allSubj),2), : ) = [];

% Save pruned output into xslx file
filename = 'output.xlsx';

% Total Datatable
xlswrite(filename,allSubj,pruned_output,'A2');
xlswrite(filename,headers,pruned_output,'A1');

% Subject Datatables
for iSubj = 1:size(outGrp,1)
    xlswrite(filename,outGrp{iSubj,1},iSubj+pruned_output,'A2');
    xlswrite(filename,headers,iSubj+pruned_output,'A1');
end

% From this point, the 'outGrp' variable represents the pruned output on
% which all total and individual subject calculations will be conducted.

%% PrePost Exercise Comparisons
% After pruning unusable studies, the remaining studies are matched with 
% their corresponding pre- or post-exercise study. Post-exe/Pre-exe ratios
% are then calculated and output in a static location as an additional set
% of columns on the individual subject sheets. 
%
% *Dependent on variable order, requires paired studies
%
% *Requires pre-occlusion plateau, ratio of max perfusion to pre-occlusive 
% plateau, and average temperature columns to correspond with their
% variables

fprintf('Completing Post / Pre Comparisons...\n');

% If a subject is removed find unique initials again
% Rebuild modified grpData for unique function call
lengthNew = zeros(1,1);
for iSubj = 1:size(outGrp,1);
    lengthNew = lengthNew + size(outGrp{iSubj,1},1);
end
grpDataMod = cell(0,numvars);
for iSubj = 1:size(outGrp,1);
    grpDataMod = vertcat(grpDataMod,outGrp{iSubj,1}); %#ok<AGROW>
end
[uniqueInit,~,~] = unique(grpDataMod(:,1),'stable');

% Save output into structure file for pre/post calculations
for iSubj = 1:size(outGrp,1)
    % Locate where pre studies start and where post studies start
    [uniquePP,ppStart,~] = unique(outGrp{iSubj,1}(:,2),'stable');
    % Loop for N-1 groups
    for iGrp = 1:length(ppStart)-1
        startGrp = ppStart(iGrp);
        endGrp   = ppStart(iGrp+1);
        currGrp  = outGrp{iSubj,1}(startGrp:endGrp-1,:);
        currInit = char(uniqueInit(iSubj,1));
        currPP   = char(outGrp{iSubj,1}(ppStart(iGrp),2));
        OutGrpPP.(currInit).(currPP) = currGrp;
    end
    % Last Group
    startGrp = ppStart(end);
    endGrp   = size(outGrp{iSubj,1},1);
    currGrp  = outGrp{iSubj,1}(startGrp:endGrp,:);
    currInit = char(uniqueInit(iSubj,1));
    currPP   = char(outGrp{iSubj,1}(ppStart(end),2));
    OutGrpPP.(currInit).(currPP) = currGrp;
end

% PrePost vars and headers; vars must be specified
ppVars = [colPlat1, colPMaxOvrPlat1, colAvgTemp];
ppHeaders = {'Post/Pre-plat1','Post/Pre-pMax/plat1','Post/Pre-avgTemp'};

structInit = fieldnames(OutGrpPP);
for iSubj = 1:size(outGrp,1); % Loop for ea. subject    
    currInit = char(structInit(iSubj));
    for iStudy = 1:size((OutGrpPP.(currInit).Post),1) % Loop for ea. study
        for iVar = 1:size(ppVars,2) % Loop for ea. variable
            % Select Data & Convert to double datatype
            ppGrp1 = cell2mat(OutGrpPP.(currInit).Post(iStudy,ppVars(iVar)));
            ppGrp2 = cell2mat(OutGrpPP.(currInit).Pre(iStudy,ppVars(iVar)));
            ppCalc = ppGrp1/ppGrp2;
            % Save
            OutGrpPP.(currInit).PP(iStudy,iVar) = ppCalc;
        end
    end
end

% To place the PP comparison data block 1 column away from the
% LDF_smooth processing data block, convert the number of columns (aka the
% number of processed variables) into the excel alphabetical column format.
% This is required because the 'xlswrite' function does not make column
% number conversions (e.g. column # 32 == AH).

column = numvars +2; % skip a column after last data block value

dict = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',...
    'p','q','r','s','t','u','v','w','x','y','z'];

dictSize = size(dict,2);

if column > size(dict,2)
    char1    = upper(dict(floor(column/dictSize)));
    char2    = upper(dict(rem(column,dictSize)));
    blockPos = strcat(char1,char2);

elseif column == dictSize
    blockPos = upper(dict(floor(column/dictSize)));

else
    char1 = upper(dict(rem(column,dictSize)));
    blockPos = char1;

end

% Generate Output
for iSubj = 1:size(outGrp,1);
    currInit = char(structInit(iSubj));
    currPP = OutGrpPP.(currInit).PP;
    xlswrite(filename,currPP,iSubj+pruned_output,strcat(blockPos,'2'));
    xlswrite(filename,ppHeaders,iSubj+pruned_output,strcat(blockPos,'1'));
end

%% Dynamically name output.xslx sheets via ActiveX control
% After saving to output.xlsx is complete, label the sheets in output.xlsx
% to match subject initials and 'total' data tables.
structInit = fieldnames(OutGrpPP);
e = actxserver('Excel.Application'); % open ActiveX server
ewb = e.Workbooks.Open([pwd '\output.xlsx']);

% Summary Sheets Loop
summarySheets = {'genTime','raw','prunedTotal'};
for iSummSheets = 1:length(summarySheets)
    ewb.Worksheets.Item(iSummSheets).Name = summarySheets{iSummSheets};
end

% Subject Loop
for iSubj = 1:size(outGrp,1)
    ewb.Worksheets.Item(iSubj+pruned_output).Name = char(structInit(iSubj));
    ewb.Save
end
ewb.Close(false);
e.Quit

%% Gen_Perm calculations and comparisons
% Generate the filtered mean analysis using pre- & post-exercise-matched 
% study data. Due to the amount of new information generated, the absolute 
% percentage difference of the filtered mean is placed into a new 
% workbook (called 'outputGenPerm.xlsx'). The other associated data, such 
% as the absolute differences and the filtered means themselves, are 
% saved into a structure variable in the Matlab workspace (called 
% 'OutGrpPP'). An additional for-loop around the subject loop could be set
% up to generate additional workbooks with this data.
%
% *Dependent in variable order and requires removal of extraneous
% studies
%
% *Requires time from end-occlusion to max perfusion, max perfusion,
% perfusion velocity, ratio of max perfusion over pre-occlusive plateau,
% and pre-occlusive pleateau columns to correspond with their
% variables
fprintf('Completing filtered means analysis...\n');

sampleSize = 1000;
numStudiesLimit = 3;

gPermVars = [colEndOccToMax,colPMax,colPVel,colPMaxOvrPlat1, colPlat1];
gPermHeaders = {'EndOccToMax','pMax','pVel','MaxOverPlat1','plat1'};

structInit = fieldnames(OutGrpPP);
structState = {'Post','Pre'};

% Pre-Allocate
genPermOutArr = cell(size(fieldnames(OutGrpPP),1)-1,1);
numGPermVars = size(gPermVars,2);

maxBlockWidth = 5; % guesstimate for preallocation
maxSummStatHeight = 11; % guesstimate for preallocation

for iSubj = 1:length(genPermOutArr)
    genPermOutArr{iSubj,1} = cell(sampleSize*2+ maxSummStatHeight, ...
                                numGPermVars*maxBlockWidth+numGPermVars-1);
end % Pre-allocation currently not 100% accurate and extensible

% Generate filtered means analysis data
for iSubj = 1:size(fieldnames(OutGrpPP),1)-1 % loop for ea. subject
    currInit = char(structInit(iSubj));
    
    % Loop for ea. variabale
    for iVar = 1:size(gPermVars,2); 
        currGPermVar = gPermVars(iVar);
        
        % Loop for ea. pre-exe or post-exe state
        for iPreOrPost = 1:size(structState,2) 
            currState = char(structState(iPreOrPost));
            currHeader = cell2mat(gPermHeaders(iVar));
            data = cell2mat(OutGrpPP.(currInit).(currState)(:,currGPermVar));
            
            % Skip subjs where # of studies are insuff
            if size(data,1) <=2 
                continue
            else
            
            % Call filtMeansSimulation script for AbsPerDiff Calculations%
            filtMeansSimulation
            genPermArr.randsmplOut = genPermArr.AbsPerDiff;
            genPermArr.randMean    = genPermArr.AbsNumStudyMeans;
            genPermArr.randStDev   = genPermArr.AbsStDev;
            genPermArr.randStErr   = genPermArr.AbsStErr;

            end % End conditional
            
            % Save headers
            OutGrpPP.(currInit).genPerm.(currState).(currHeader) = genPermArr;
           
            % Re-organize data into blocks for excel output
            blockHeight = size(num2cell(genPermArr.randsmplOut),1);
            blockWidth = size(num2cell(genPermArr.randsmplOut),2);
            
            % Rearrange Sample data
            randsmplOutYPos = (iPreOrPost-1)*blockHeight+(iPreOrPost);
            randsmplOutXPos = (iVar-1)*blockWidth+(iVar-1);
            genPermOutArr{iSubj,1}( ...
                (1+randsmplOutYPos): ... % y pos for data block
                (blockHeight+randsmplOutYPos), ...
                (1+randsmplOutXPos): ... % x pos for data block
                (blockWidth+randsmplOutXPos) ...
                ) = num2cell(genPermArr.randsmplOut);
            
            % Capture generated summary statistics
            gPermSummaryStat = {genPermArr.randMean,...
                                genPermArr.randStDev,...
                                genPermArr.randStErr};
             
            % Loop for ea. extra pre and post data block
            for iSummStat = 1:length(gPermSummaryStat)
                
                % Starting after randPerm data, skip a line after ea block
                blockYPos = blockHeight*2+3+iPreOrPost+(iSummStat-1)*3;
                blockXPos = randsmplOutXPos;
                genPermOutArr{iSubj,1}( ...
                    blockYPos:blockYPos, ...
                    (1+blockXPos):(blockWidth+blockXPos) ...
                    ) = num2cell(gPermSummaryStat{iSummStat});
            end % end pre/post data block loop
        end % end exe-state loop
        
        % Variable Labels
        headerYPos = 1;
        headerXPos = randsmplOutXPos;
        genPermOutArr{iSubj,1}(1,1+headerXPos) = gPermHeaders(iVar);    
        
        % TTEST (unpaired, homoscedastic)
        for iCol = 1:blockWidth
            % Select Column
            grp1 = genPermOutArr{iSubj,1}(...
                (1+1):(blockHeight+1),(iCol));
            grp2 = genPermOutArr{iSubj,1}(...
                (1+(blockHeight+2)):((blockHeight+2)+blockHeight),(iCol));
            
            % Convert to double
            grp1 = cell2mat(grp1);
            grp2 = cell2mat(grp2);
            
            % 2-sample T-Test (paired-sample, homoscedastic, 2-tailed)
            [h,p] = ttest(grp1,grp2);
            
            % Write result to array
            genPermOutArr{iSubj,1}...
                (blockYPos+2,(1+randsmplOutXPos+(iCol-1))) = {p};
        end    
    end
end


% Labels for excel sheet (did not use recursive fcn for-loop for clarity)
gpLabels = cell(63,1);
n = 1; % Start
sS = sampleSize+1; % Interval
gpLabels(n) = {'AbsPerDiff Vars'};
gpLabels(n+3) = {'Post'};
gpLabels(n+3+sS) = {'Pre'};
gpLabels(n+3+(sS*2)) = {'PostMean'};
gpLabels(n+3+(sS*2)+1) = {'PreMean'};
gpLabels(n+3+(sS*2)+3) = {'PostStDev'};
gpLabels(n+3+(sS*2)+4) = {'PreStDev'};
gpLabels(n+3+(sS*2)+6) = {'PostStErr'};
gpLabels(n+3+(sS*2)+7) = {'PreStErr'};
gpLabels(n+3+(sS*2)+9) = {'T-Test'};

fprintf('Saving filtered means analysis...\n');
filename2 = 'outputGenPerm.xlsx';

% Export current time
currTime = datestr(now);
currTime = {'Gen on'; currTime};
xlswrite(filename2,currTime,1,'A1');

% Todo: combine labels and tables for less xlswrite calls

% Write output to outputGenPerm.xslx
for iSubj = 1:size(fieldnames(OutGrpPP),1)-1;
    currGenPerm = genPermOutArr{iSubj,1};
    xlswrite(filename2,currGenPerm,iSubj+1,'B3');
    xlswrite(filename2,gpLabels,iSubj+1,'A1');
end

% Dynamically name outputGenPerm.xslx sheets via ActiveX control
structInit = fieldnames(OutGrpPP);
e = actxserver('Excel.Application'); % open ActiveX server
ewb = e.Workbooks.Open([pwd '\outputGenPerm.xlsx']);
ewb.Worksheets.Item(1).Name = 'genTime';

% Subject Loop
for iSubj = 1:size(fieldnames(OutGrpPP),1)-1
    ewb.Worksheets.Item(iSubj+1).Name = char(structInit(iSubj));
    ewb.Save
end 
ewb.Close(false);
e.Quit

% save('output.mat','outputArray');

%% Resampling Analysis
% Not implemented


%% Summary statistics
% After most of the data is generated, generate the mean, st.deviation,
% and standard error for all subjects, individual subjects, and all pre- or
% post-exercise studies for individual subjects. A T-test(paired, 2-tailed, 
% homoscedastic) is run to compare pre- and post-exercise studies in 
% individual subjects.
%
% To do this, we first set up a double data-type table for the summary 
% statistics. We can then assign the results into a 'SummaryStats' field in 
% the 'OutGrpPP' structure and write the results into output.xlsx.
% 

fprintf('Saving summary statistics...\n');

% All subjects
% Pre-allocate Data tables/Convert data table to double for number functions
OutGrpPP.SummaryStats.data = grpDataMod;
OutGrpPP.SummaryStats.dataSumm(:,4:numvars) = cell2mat(grpDataMod(:,4:numvars));
summStatsTotal = zeros(3,numvars);

% Flag total, pre-exe, post-exe studies
flagAllStudies  = ((OutGrpPP.SummaryStats.dataSumm(:,4)>0));
flagPostStudies = (strcmp(grpDataMod(:,2),'Post'));
flagPreStudies  = (strcmp(grpDataMod(:,2),'Pre'));
summStatsSet    = {flagAllStudies,flagPostStudies,flagPreStudies};

% Assign summary stats
statblock = zeros(4*size(summStatsSet,2),numvars);
for x = 1:size(summStatsSet,2)
    exeStatBlock = zeros(3,numvars); % 3 summary stats
    for column = 4:size(grpDataMod,2)
        data = OutGrpPP.SummaryStats.dataSumm(summStatsSet{x},column);
        exeStatBlock(1,column-3) = mean(data);
        exeStatBlock(2,column-3) = std(data);
        exeStatBlock(3,column-3) = std(data)/sqrt(numel(data));
    end
    
    % Arrange stat blocks: Total -> Post-exe -> Pre-exe
    indBlock = x+(x-1)*3;
    statblock(indBlock:indBlock+2,:) = exeStatBlock;
end

% T-test
if x == 3
    
    % Initialize block
    tTestBlock = zeros(1,numvars); 
    
    % Loop for input into column
    for column = 4:size(grpDataMod,2)
        data1 = OutGrpPP.SummaryStats.dataSumm(summStatsSet{2},column);
        data2 = OutGrpPP.SummaryStats.dataSumm(summStatsSet{3},column);
        
        % paired-sample, homoscedastic, 2-tailed)
        [~,tTestBlock(1,column-3)] = ttest(data1,data2);
    end % end column loop
    
    % Merge into summary statblock
    indTTestBlock = x + (x - 1) * 3 + 4;
    statblock(indTTestBlock,:) = tTestBlock;
end

statblock(statblock == 0) = NaN;

% Subject sheets; requires 'outGrp' variable
% Pre-allocate Pre & Post tables
structInit = fieldnames(OutGrpPP);
structState = {'Post','Pre'};

summStatsPP = cell(size(outGrp,1),1);
for iSubj = 1:size(outGrp,1)
% 3 vars for total, pre- and post-exe + TTest = 3*3+1+ 3 lines skipped = 13
summStatsPP{iSubj,1} = zeros(13,numvars); 
end

for iSubj = 1:size(outGrp,1) % uses outGrp size which ignores SummaryStats
    currInit = char(structInit(iSubj));
    % Individual Totals
        % Set up data table for total summary stats
        blockSizeX = size(outGrp{iSubj,1},2);
        blockSizeY = size(outGrp{iSubj,1},1);
        OutGrpPP.(currInit).SummaryStats.Total.data = outGrp{iSubj,1};
        OutGrpPP.(currInit).SummaryStats.Total.dataSumm(1:blockSizeY,4:blockSizeX) = ...
          cell2mat(outGrp{iSubj,1}(1:blockSizeY,4:blockSizeX));
      
        % Assign summary stats
        for b = 4:size(OutGrpPP.(currInit).SummaryStats.Total.dataSumm,2);
            
            data = OutGrpPP.(currInit).SummaryStats.Total.dataSumm(:,b);
            currMean  = mean(data);
            currStd   = std(data);
            currStErr = currStd/sqrt(size(data,1));

            summStatsPP{iSubj,1}(1,b) = currMean;
            summStatsPP{iSubj,1}(2,b) = currStd;
            summStatsPP{iSubj,1}(3,b) = currStErr;
                    
            OutGrpPP.(currInit).SummaryStats.Total.mean(b) = currMean;
            OutGrpPP.(currInit).SummaryStats.Total.stDev(b) = currStd;
            OutGrpPP.(currInit).SummaryStats.Total.stErr(b) = currStErr;
        end % end column loop
        
    % Individual Pre- and Post-exercise
    for iPreOrPost = 1:size(structState,2)
        currState = char(structState(iPreOrPost));
        
        % Set up data table for summary stats
            OutGrpPP.(currInit).SummaryStats.(currState).dataSumm(:,4:numvars) = ...
                cell2mat(OutGrpPP.(currInit).(currState)(:,4:numvars));
        
        % Assign summary stats
        for iCol = 4:size(OutGrpPP.(currInit).SummaryStats.(currState).dataSumm,2)
            
            % Select Data
            data = OutGrpPP.(currInit).SummaryStats.(currState).dataSumm(:,iCol);
            
            % Number functions
            currMean  = mean(data);
            currStd   = std(data);
            currStErr = currStd/sqrt(size(data,1));
            
            % Arrange stat block: Skip a line after individual totals(+4), 
            % then skip a line for pre-exercise studies ((iPreOrPost-1)*4)
            row = (iPreOrPost-1)*4+4;
            
            % Mean -> StDev -> StErr
            summStatsPP{iSubj,1}(1+row,iCol) = currMean;
            summStatsPP{iSubj,1}(2+row,iCol) = currStd;
            summStatsPP{iSubj,1}(3+row,iCol) = currStErr;
            
            OutGrpPP.(currInit).SummaryStats.(currState).mean(iCol)  = currMean;
            OutGrpPP.(currInit).SummaryStats.(currState).stDev(iCol) = currStd;
            OutGrpPP.(currInit).SummaryStats.(currState).stErr(iCol) = currStErr;
        end % end column loop
    end % end pre-exe/post-exe loop
    
    % Assign T-test data
    for iCol = 4:size(OutGrpPP.(currInit).SummaryStats.(currState).dataSumm,2)
        
        % Select columns
        grp1 = OutGrpPP.(currInit).SummaryStats.Post.dataSumm(:,iCol);
        grp2 = OutGrpPP.(currInit).SummaryStats.Pre.dataSumm(:,iCol);
        
        % 2-sample T-Test (paired-sample, homoscedastic, 2-tailed)
        [h,p] = ttest(grp1,grp2);
        
        % Write result to structure
        OutGrpPP.(currInit).SummaryStats.PrePostTTest(iCol) = {p};
        
        % Write result to summStats table shifted down 2 cells from
        % previously generated summary stats
        summStatsPP{iSubj,1}(3+row+2,iCol) = p;
    end % end column loop
end % end subject summ stats loop


% Remove zeros from summary stat cells
summStatsTotal(summStatsTotal == 0) = NaN;
for iSubj = 1:size(outGrp,1)
  summStatsPP{iSubj,1}(summStatsPP{iSubj,1}==0) = NaN;
end

% Labels for excel sheet
summStatsLabels   = cell(13,1);
summStatsHeaders1 = {'Total ', 'Post ', 'Pre '};
summStatsHeaders2 = {'Mean', 'St.Dev', 'St.Err'};

% Offset labels
for a = 1:size(summStatsHeaders1,2)
    currHeader1 = summStatsHeaders1(a);
    
    for b = 1:size(summStatsHeaders2,2);
        currHeader2 = summStatsHeaders2(b);
        combHeader  = strcat(currHeader1,currHeader2);
        currOffset  = 1+(a-1)*4;
        summStatsLabels(currOffset+b-1,1) = cellstr(combHeader);
    end
end
summStatsLabels(13) = {'T-Test'};

% Todo: combine labels and tables for less xlswrite calls

% Excel output
row = num2str(size(allSubj,1)+3);
xlswrite(filename,num2cell(statblock),pruned_output,strcat('D',row))
xlswrite(filename,summStatsLabels,pruned_output,strcat('A',row))

for iSubj = 1:size(outGrp,1)
    row = num2str(size(outGrp{iSubj,1},1)+3); 
    
    % skip a line after data block (+2), shift down 1 row for headers (+1)
    xlswrite(filename,summStatsPP{iSubj,1},iSubj+pruned_output,strcat('A',row));
    xlswrite(filename,summStatsLabels,iSubj+pruned_output,strcat('A',row));
end

%% Save OutGrpPP
fprintf('Saving structure variable...\n');
filename3 = 'MatlabStructureVarOutput.mat';

% Export current time
currTime = datestr(now);

% Make folder + set current Dir if they don't exist
currDir = pwd; 
if exist([pwd '\Structure Var Output'],'dir') == 0
    mkdir('Structure Var Output')
end

% Save into new folder
cd([currDir '\Structure Var Output'])
save(filename3, 'OutGrpPP','currTime','matfiles','outputArray','numvars'...
    ,'outputEvalVars','inputOccTimings','allSubj');
cd(currDir)

%% Generate Raw tracings
% For comparisons with smoothed data analysis.
fprintf('Generating raw tracings for reference...\n');
generate_raw_tracings

%% Cleanup

% Remove all variables directly unrelated to viewing the final output and
% processing errors.
clearvars -except grpData...
                  outGrp...
                  OutGrpPP...
                  matfiles...
                  in_time_restr...
                  inStartOcc...
                  inEndOcc...
                  in_minpkdis...
                  in_perf_limit...
                  err...
                  outputArray...
                  outputEvalVars...
                  numvars...
                  inputOccTimings
disp('done!')

% Turn off error Logging
diary off