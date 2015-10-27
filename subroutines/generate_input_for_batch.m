% Generates an excel sheet for collecting occlusion data manually
% Remember to import input.xlsx before starting the batch script
matfiles = dir('*.mat');
matarr = struct2cell(matfiles);
matnames = transpose(matarr(1,:));
headers = {'names' 'inStartOcc' 'inEndOcc'};

% Structure File Output
LDF_input = struct(...
    'names',matnames,'inStartOcc',[],'inEndOcc',[]);

% Excel Output
filename = 'input.xlsx';
xlswrite(filename,headers,2,'A1');
xlswrite(filename,matnames,2,'A2'); %new directory on 2nd sheet

e = actxserver('Excel.Application'); % open ActiveX server
ewb = e.Workbooks.Open([pwd '\' filename]);
ewb.Worksheets.Item(1).Name = 'input';
ewb.Worksheets.Item(2).Name = 'new ordered matfiles';
ewb.Save
ewb.Close(false);
e.Quit

clearvars matarr matnames headers filename sheet trans xlRange ewb e

fprintf('Generated input.xlsx in current folder.\nNew data is on sheet 2.\n')