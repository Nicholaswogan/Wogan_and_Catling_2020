% This program searches through the cell array 'databaseB' for name matches
% with the parameter string 'S in the first index. 

function index = searchDataB(S)
    global databaseB;                        % Load the database information
    list = regexp(databaseB(:,1),S,'once','freespacing');  % Create an array with the number
                                            % of occurences in it.
    index = -1;
      
    for x=1:length(list)                    % Check each of the occurences
        if length(list{x,1})>0 && (list{x,1}(1)==1)  % If the number of them that start with the search term is non-zero 
            index=x;                        % Remember the last one.
            break
        end
    end
    clear('list','x','hits')                % Clear all non-essential variables.
end