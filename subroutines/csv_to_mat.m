%Convert .csv doppler studies in current directory to .mat
csvfiles = dir('*.csv'); 
for i = 1:length(csvfiles)
    study = mat2str(csvfiles(i).name);
    b = importdata(csvfiles(i).name);
    time =(b.data(:,1)); %#ok<*NASGU> % Variables will be used elsewhere
    temp = (b.data(:,2));
    pressure = (b.data(:,3));
    perfusion = (b.data(:,4));
    study = strrep(study,'.csv','');
    study = strrep(study,'''','');
    studyNameOutput = strsplit(study);
    initial = cell2mat(studyNameOutput(1));
    prePost = cell2mat(studyNameOutput(2));
    date = cell2mat(studyNameOutput(3));
    save(study,'study','initial','prePost','date','time','temp','pressure','perfusion');
    clear study b time temp pressure perfusion study studyNameOutput initial prePost date
end
