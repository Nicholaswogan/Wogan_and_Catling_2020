% This program loads the aqueous phase database in the file 'F' into the cell array
% 'databaseC'

F = 'sprons96_edited2.dat';      % same as original except CO renamed to find it more easily
data = [];              % Variable used for holding each species' data temporarily
n = 0;                  % Current number of species in the database.
fid = fopen(F);         % Opens the data file and stores the I/O id number
global databaseC;
databaseC = cell(1,2);

trash = '';

while isempty(strfind(trash, '      aqueous species'))           % Removes all the intro comments
    trash = fgetl(fid);
end
trash = fgetl(fid);


line = fgetl(fid);      % Check the first line.

while ~isempty(line)      % While there are more lines to read. . .
     
    n = n + 1;          % Increase the number of species in the database.
    
    databaseC(n,1) = {strtrim(line(20:end))};    % Place the species name in the first column of the database
    for temp = 1:3
        line = fgetl(fid);                          % Get the next line.
    end
    
    coefA = textscan(line(5:end),'%12f',3);        % Store the range of temperature values that this data applies to
    line = fgetl(fid);                      % Get the next line
    coefB = textscan(line(2:end),'%12f',4);        % Store the coefficients of the first line of data
    line = fgetl(fid);                      % Get the next line
    coefC = textscan(line(2:end),'%12f',3);        % Store the coefficients of the second line of data
    coefD = textscan(line(39:end),'%2f',1);
    
    data = [coefA{1};coefB{1};coefC{1};coefD{1}]; % Stack the three sets of data obtained into one column and store it into the data matrix.


    databaseC(n,2) = {data};                     % Add the matrix to the database location across from the species name.
    line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
end

databaseC(n,1)
data

save databaseC
clear('F','data','fid','x','trash','line','temperatures','range','coefA','coefB','coefC','coefD','temp','n');
% Clear all variables associated with the program except the database and
% the total number of species in the database.