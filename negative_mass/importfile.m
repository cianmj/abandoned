function importfile(fileToRead1)
% Import the file
rawData1 = importdata(fileToRead1);

% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.
[unused,name] = fileparts(fileToRead1);
newData1.(genvarname(name)) = rawData1;

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

