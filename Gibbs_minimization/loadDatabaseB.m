% This program loads the SECOND database in the file 'F' into the matrix 'databaseB'.

F = 'fugacityCoefficientVariables.txt';      % The name of the file holding the data
fid = fopen(F);         % Opens the data file and stores the I/O id number
global databaseB
databaseB = cell(0,4);

line = fgetl(fid);      % Delete the Trash Line
line = fgetl(fid);      % Check the first line.

while ischar(line)      % While there are more lines to read. . .
    
    while length(line)<5    % Discard any lines without any information    
        line=fgetl(fid);    % ==0 would also work, but this handled      
    end                     % an odd space which stopped the program
    
    line = regexprep(line,line(1),'');                    % Removes an unidentified character from the string
    databaseB = [databaseB ;regexp(line,'\t','split')];   % Place the species name in the first column of the database
    line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
end
clear('F','fid','x','trash','line');
% Clear all variables associated with the program except the database and
% the total number of species in the database.