%Convert .csv doppler studies in current directory to .mat
csvfiles = dir('*.csv'); 

fprintf('Converting  .csv to .mat...\n');

for i = 1:length(csvfiles)
    
    % Progress
    fprintf(' %s \n',csvfiles(i).name);
    
    % Import csv data
    study = mat2str(csvfiles(i).name);
    b = importdata(csvfiles(i).name);
    time =(b.data(:,1)); %#ok<*NASGU> % Variables will be used elsewhere
    temp = (b.data(:,2));
    pressure = (b.data(:,3));
    perfusion = (b.data(:,4));
    
    % Parse csv filename into study identifier variables
    study = strrep(study,'.csv','');
    study = strrep(study,'''','');
    studyNameOutput = strsplit(study);
<<<<<<< HEAD
    
    % Study Identifier variables
    patn    = cell2mat(studyNameOutput(1));
    initial = cell2mat(studyNameOutput(2));
    visit   = cell2mat(studyNameOutput(3));
    prePost = cell2mat(studyNameOutput(4));
    date    = cell2mat(studyNameOutput(5));
    
    % Save mat file
    save(study,'study','patn','initial','visit','prePost','date','time','temp','pressure','perfusion');
    
    % Cleanup
    clear study b time temp pressure perfusion study studyNameOutput initial prePost date patn visit
=======
    initial = cell2mat(studyNameOutput(2));
    prePost = cell2mat(studyNameOutput(4));
    date = cell2mat(studyNameOutput(5));
    save(study,'study','initial','prePost','date','time','temp','pressure','perfusion');
    clear study b time temp pressure perfusion study studyNameOutput initial prePost date
>>>>>>> origin/master
end

fprintf('done! \n');
