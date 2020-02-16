% This program searches through the cell array 'databaseD' for name matches
% with the parameter string 'S in the first index.

function [index1,index2] = searchDataD(S)
    global databaseD;                        % Load the database information

    for x=1:length(databaseD(:,1))                    % Check each of the occurences
        if strcmp(S,databaseD(x,1))==1  % If the number of them that start with the search term is non-zero 
            index1=str2double(databaseD(x,2));                        % Remember the last one.
            index2=str2double(databaseD(x,3));
            break
        end
    end
end