% This program searches through the cell array 'databaseC' for name matches
% with the parameter string 'S in the first index. If it finds more than
% one option, it returns -2. If it doesn't find any options, it returns -1.
% Otherwise, it returns the index of the species with name S.

function index = searchDataC(S)
    global databaseC;                        % Load the database information
    list = regexp(databaseC(:,1),S,'once','freespacing');  % Create an array with the number
                                            % of occurences in it.
    index = -1;
    temp = 0;  
    for x=1:length(list)                    % Check each of the occurences
        if ~isempty(list{x,1}) && (list{x,1}(1)==1)  % If the number of them that start with the search term is non-zero 
            index=x;                        % Remember the last one.
            temp = temp+1;
        end
    end
    if temp==0                      % If there are no indexes
        index=-1;                           % set the return value to -1
    end
    if temp>1                       % If there is more than one index
        index=-2;                           % set the return value to -2.
    end
    
    clear('list','x','temp')                % Clear all non-essential variables.
end