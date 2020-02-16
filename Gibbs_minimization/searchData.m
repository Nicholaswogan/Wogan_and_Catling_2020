% This program searches through the cell array 'database' for name matches
% with the parameter string 'S in the first index. If it finds more than
% one option, it returns -2. If it doesn't find any options, it returns -1.
% Otherwise, it returns the index of the species with name S.

function index = searchData(S)
    global database;                        % Load the database information
    list = regexp(database(:,1),S,'once');  % Create an array with the number
                                            % of occurences in it.                                        
    hits = [];                              % Array of all successful matches
    for x=1:length(list)                    % Check each of the occurences
        if length(list{x,1})>0 && (list{x,1}(1)==1)  % If the number of them that start with the search term is non-zero 
            index=x;                        % Remember the last one.
            hits = [hits index];            % Add this to the list of successful finds
        end
    end
    
    if length(hits)==0                      % If there are no indexes
        index=-1;                           % set the return value to -1
    end
    if length(hits)>1                       % If there is more than one index
        index=-2;                           % set the return value to -2.
    end
    clear('list','x','hits')                % Clear all non-essential variables.
end