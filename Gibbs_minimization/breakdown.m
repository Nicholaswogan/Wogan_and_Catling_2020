% Returns a cell array of all of the elements in a given species in the first row and the number of each element in the second row.

function n=breakdown(str)
    %str:            The species being parsed (eg: H2O, HCl, etc.)
    n = cell(0);%    The cell array Elements and the number of each element in the species being parsed.
    
    %                                   This code strips off any non-alphanumeric characters in the species name
    temp = isstrprop(str,'alphanum');%  Returns an array of ones and zeros. Anything that is not a number of letter is a zero.
    for x=1:length(str)%                This loop finds the first non-alphanumeric character, then cuts the str at that point so that
        if temp(x)==0%                       the original string only contains alphanumeric characters (this is useful for certain species
            str = str(1:x-1);%               like OH(-)   )
            break
        end
    end
    
    %Once the str is only alphanumeric, the parser seperates it.
    n=parser(str);
        
end

    % The parser returns a cell array as described above. NOTE: The function is recursive, meaning that it calls itself.
function c=parser(str)

    % In the rare case that the passed str is null, return a blank cell array.
    if isempty(str)
        c=cell(0);
    else                   % Otherwise, remember the first character and everything else.
        front = str(1);
        rest = str(2:end);
        if isstrprop(front,'upper')                          % If the first character is NOT an uppercase letter, the species was improperly named.
            if ~isempty(rest) && isstrprop(rest(1),'lower')  % If the passed str is more than one character and the second character is
                front = [front rest(1)];                     % lowercase, attach it to the previous letter (C l becomes Cl)
                rest = rest(2:end);                          % Also update the 'rest'
            end
            if isempty(rest)                                 % If there are no more characters, return the cell array {element name stored
                c={front;1};                                 % in 'front' ; 1}
            elseif isstrprop(rest(1),'digit')                % If there ARE more characters, and the next char is a number, figure out the
                d=grabNumbers(0,rest);                       % number, and then make 'rest' become everything after the number.
                rest = rest(1+ceil(log10(d+1)):end);         % log10(d + 1) finds the number of characters in the number of elements.
                if isempty(rest)                             % Now if there are no more characters, return the cell array {element; number}
                    c={front;d};
                else                                         % Otherwise, parse everything else
                    c=parser(rest);
                    place=findElement(front,c);              % If the current element is one of the 'everything else' don't add a new
                    if place==-1                             % column, just increase the number of that element (this handles CH2OH)
                        c=[{front;d} c];                     % Otherwise, just add the current element to 'everything else'
                    else
                        c{2,place}=(c{2,place}+d);           
                    end
                end
            else
                c=parser(rest);                              % Since there are more characters, but they are different elements, parse them.
                place=findElement(front,c);                  %  Then, add that element to the cell array. If the element is not already in 
                if place==-1                                 %  the cell array, then add it like normal to the cell array. Otherwise,
                    c=[{front;1} c];                         %  find the location of the other element and add it to that location.
                else
                    c{2,place}=(c{2,place}+1);
                end 
            end
        else
            disp(['Illegal Name: ' str])                     % Since the species was improperly named, throw an error.
            throw(MException('Catling:InvalidSpecies','Invalid species passed to breakdown'))
        end
    end
end

% This function returns the number at the front of a string (used to find the number of elements). NOTE: This function is recursive.
function d=grabNumbers(val,str)
    if isempty(str)                    % If there are no more numbers, then there can't be any more numbers to pick up.
        d=val;                         % In that case, just return the answer.
    elseif isstrprop(str(1),'digit')   % Otherwise, if there are more characters, and those characters are number, add them.
        d=grabNumbers(str2num(str(1))+10*val,str(2:end));
    else
        d=val;                         % If there are more characters, but they are letters, then the function is done. Just return the answer.
    end
end