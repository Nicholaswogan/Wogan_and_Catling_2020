% Returns the location in the passed horizontal string cell array 'a' of the passed string 'str'. Returns -1 if not in the cell array.

function n=findElement(str,cellarray)
    % str: A string indicating which element is being searched for in the cell array
    % cellarray: A cell array containing elements and the number of each element in a given species (this cell array is typically the type
    %               that would be returned from breakdown
    n = -1;                              % The location of the element in the given cell array.
    [x,y] = size(cellarray);
    for i=1:y                            % For each element in the given cell array.
        if strcmp(str,cellarray{1,i})    % If the string in the first column is the given element (ie, the element is already in the cell array)
            n=i;                         % then set the location of that string to be n and stop executing the program.
            break;
        end
    end
    clear('x','y') % Clear variables from memory.
end
        
