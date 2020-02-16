%Takes a cell array of species names and creates the coefficient matrix for it.
%NOTE: The top row has all the elements in it.
% Precondition: names must have at least one species.

function a=makeAs(names)
    if ~iscellstr(names)           % If only one species is being parsed for the coefficient matrix, the answer will just be the result
        a=breakdown(names);        %     from breakdown.
    else
        a=breakdown(names{1});     % If there are more than one species, then start by running breakdown on the first one.
        names = names(2:end);      % Use this as the starting point for the rest. Remove the first thing from the things that still
                                   % need to be parsed.
        for i = 1:length(names)    % For each species still needing to be parsed.
            [x,y] = size(a);       %        find the size of a.
            a = [a ; num2cell(zeros(1,y))];% Then add an additional row to a (for the new species)
            new = breakdown(names{i});     % Now run breakdown on that species.
            [x2,y2] = size(new);           % Find the size of the result.
            for j = 1:y2                   % For each element in the result.
                f = findElement(new{1,j},a);% See if the element is already in the coefficient matrix.
                if f==-1                    % If it is not, add a new column to the coefficient matrix and stick that element on the end.
                    a = [a [new{1,j} ; num2cell(zeros(x-1,1));new{2,j}]]; % Filling in the gap for the rest of the species (the other rows)
                else                        % If it is already in the list, just add the element to its proper column. 
                    a{end,f}=new{2,j};
                end
            end
        end
    end
    
end