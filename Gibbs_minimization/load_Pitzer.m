% This script loads the database of Pitzer binary interaction parameters
% from their database file and saves them in the cells B0, B1, B2, and C0.

F = 'pitzer_coef.txt';      % The name of the file holding the data
fid = fopen(F);         % Opens the data file and stores the I/O id number
global B0 B1 B2 C0
%%% Create empty cells for B0, B1, B2, and C0:
B0 = cell(0,6);
B1 = cell(0,6);
B2 = cell(0,6);
C0 = cell(0,6);

line = fgetl(fid);      % Check the first line.

% define indices for sequentially filling out each cell
in_B0=1;
in_B1=0;
in_B2=0;
in_C0=0;
while ischar(line)      % While there are more lines to read. . .
    
    
    line = regexprep(line,line(1),' ')       % Removes an unidentified character from the string
    if strcmp(line,' B0')
        line = fgetl(fid); % initially store values for B0 cell
    elseif strcmp(line,' B1') % if you reach line defining B1 parameters, change indices to store values in B1 cell
        line = fgetl(fid);
        in_B0=0;
        in_B1=1;
    elseif strcmp(line,' B2') % if you reach line defining B2 parameters, change indices to store values in B2 cell
        line = fgetl(fid);
        in_B1=0;
        in_B2=1;
    elseif strcmp(line,' C0') % as above but for C0 cell.
        line = fgetl(fid);
        in_B2=0;
        in_C0=1;
    end

    temp_line=strsplit(line);
    temp_line=temp_line(2:end);
        
    if in_B0==1
        B0 = [B0;temp_line];   % Fill in line in B0
        line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
    elseif in_B1==1
        B1 = [B1;temp_line];   % Fill in line in B1
        line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
    elseif in_B2==1
        B2 = [B2;temp_line];   % Fill in line in B2
        line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
    elseif in_C0==1
        C0 = [C0;temp_line];   % Fill in line in C0
        line = fgetl(fid);                          % Get the next line to prepare for the while loop above.
    end
end


clear('F','fid','x','trash','line');
% Clear all variables associated with the program except the database and
% the total number of species in the database.