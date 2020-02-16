% % This program loads the 4th database in the file 'F' into the matrix 'databaseD'.

F = 'AQ_coefficients.txt';      % The name of the file holding the data
fid = fopen(F);         % Opens the data file and stores the I/O id number
global databaseD
databaseD = cell(0,3);

line = fgetl(fid);      % Delete the Trash Line
line = fgetl(fid);      % Check the first line.

while ischar(line)      % While there are more lines to read. . .
    
    regexp(line,'\t','split')
    line = regexprep(line,line(1),'');                    % Removes an unidentified character from the string
    databaseD = [databaseD ;regexp(line,'\t','split')];   % Place the species name in the first column of the database
    line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
end
clear('F','fid','x','trash','line');
% Clear all variables associated with the program except the database and
% the total number of species in the database.