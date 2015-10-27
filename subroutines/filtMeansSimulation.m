clear genPermArr

% data = rand(8,4);

if exist('sampleSize','var') == 0 && exist('numStudiesLimit','var') ==0
    numStudiesLimit = 5; % max number of studies to be averaged at once
    sampleSize      = 1000; % max number of iterations for simulation
end

% Pre-allocate
grandMean   = zeros(size(data,2),1);
means       = zeros(sampleSize,numStudiesLimit);
absPerDiffs = zeros(sampleSize,numStudiesLimit);
genPermArr  = struct('Means',[],...
                     'AbsPerDiff',[],...
                     'GrandMean',[]);
                 
numStudiesMeans = zeros(1,numStudiesLimit);
stDev = zeros(1,numStudiesLimit);
stErr = zeros(1,numStudiesLimit);

% For CV% calculations
stDevCV = zeros(1,numStudiesLimit);
stErrCV = zeros(1,numStudiesLimit);
numStudiesMeansCV = zeros(1,numStudiesLimit);

% Loop for each variable in 'data' (represented as individual columns)
for iColumn = 1:size(data,2)
    column = data(:,iColumn);
    grandMean = mean(column(:,1));
    
    % Parallel loop for increasing number of PORH tests to be averaged
    parfor iNumStudies = 1:numStudiesLimit
        
        % Loop for sample size
        for iSample = 1:sampleSize
        
        % Generate data
        currSample     = datasample(column, iNumStudies,'replace',true);
        currMean       = mean(currSample);
        currAbsPerDiff = (abs((currMean - grandMean))/...
                        grandMean)*100;
        currCV         = std(currSample)/grandMean;
                    
        % Collect data
        means(iSample,iNumStudies)       = currMean;
        absPerDiffs(iSample,iNumStudies) = currAbsPerDiff;
        CV(iSample,iNumStudies)          = currCV; % for CV% calculations
        end        
    end
    
    % Variable statistics loop for absperdiffs only
    for iNumStudies = 1:numStudiesLimit
        currAbsPerDiffs = absPerDiffs(:,iNumStudies);
        stDev(iNumStudies) = std(currAbsPerDiffs);
        stErr(iNumStudies) = std(currAbsPerDiffs)/sqrt(sampleSize);
        numStudiesMeans(iNumStudies) = mean(absPerDiffs(:,iNumStudies));
    end
        % Variable statistics loop for CV only
    for iNumStudies = 1:numStudiesLimit
        currStDev = CV(:,iNumStudies);
        stDevCV(iNumStudies) = std(currStDev);
        stErrCV(iNumStudies) = std(currStDev)/sqrt(sampleSize);
        numStudiesMeansCV(iNumStudies) = mean(CV(:,iNumStudies));
    end
    
    % Note: finding CV results in an invalid value for studies with one 
    % sample; due to SD = 0 -> 0/1 -> 0
    
    
    % Assign to structure variable
    genPermArr(iColumn).Means            = means;
    genPermArr(iColumn).AbsPerDiff       = absPerDiffs;
    genPermArr(iColumn).CV               = CV;
    genPermArr(iColumn).GrandMean        = grandMean;
    genPermArr(iColumn).AbsNumStudyMeans = numStudiesMeans;
    genPermArr(iColumn).AbsStDev         = stDev;
    genPermArr(iColumn).AbsStErr         = stErr;
    genPermArr(iColumn).CVNumStudyMeans  = numStudiesMeansCV;
    genPermArr(iColumn).CVStDev          = stDevCV;
    genPermArr(iColumn).CVStErr          = stErrCV;
    
end

