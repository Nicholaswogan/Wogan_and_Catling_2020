% This program loads the database in the file 'F' into the cell array
% 'database'.

F = 'NEWNASA.txt';      % The name of the file holding the data
data = [];              % Variable used for holding each species' data temporarily
n = 0;                  % Current number of species in the database.
fid = fopen(F);         % Opens the data file and stores the I/O id number
global database
database = cell(1,2);

for x = 1:35            % Removes all the intro comments
    trash = fgetl(fid);
end

line = fgetl(fid);      % Check the first line.

while ischar(line)      % While there are more lines to read. . .
    
    while length(line)<5    % Discard any lines without any information    
        line=fgetl(fid);    % ==0 would also work, but this handled      
    end                     % an odd space which stopped the program
    
    data = [];          % Clear the data array
    n = n + 1;          % Increase the number of species in the database.
    
    database(n,1) = textscan(line,'%16c',1);    % Place the species name in the first column of the database
    line = fgetl(fid);                          % Get the next line.
    temperatures = textscan(line,'%2n',1);      % Read the number of temperature ranges to be searched.
    
    for x = 1:temperatures{1}                   % For each temperature range . . .
        line = fgetl(fid);                      % Get the next line of data.
        range = textscan(line,'%11f',2);        % Store the range of temperature values that this data applies to
        line = fgetl(fid);                      % Get the next line
        coefA = textscan(line,'%16f',5);        % Store the coefficients of the first line of data
        line = fgetl(fid);                      % Get the next line
        coefB = textscan(line,'%16f',5);        % Store the coefficients of the second line of data
        data = [data [range{1};coefA{1};coefB{1}]]; % Stack the three sets of data obtained into one column and store it into the data matrix.
    end
                                                % Once all the matricies
                                                % are stored,
    database(n,2) = {data};                     % Add the matrix to the database location across from the species name.
    line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
end
save database
clear('F','data','fid','x','trash','line','temperatures','range','coefA','coefB','n');
% Clear all variables associated with the program except the database and
% the total number of species in the database.