% Takes a cell array of elements and the temperature in Kelvin and returns those elements in "standard form" for calculating deltaG

function cellarray=convert2standard(cellarray,T)
    [x,y]=size(cellarray);
    
    for i=1:y
        switch cellarray{1,i}                           % For each element in the cell array of elements, find the letter that they match.
            case {'Br','I','N','Cl','H','O','F'}        % If the element is diatomic, add a two afterward.
                cellarray{1,i}=[cellarray{1,i} '2 '];
            case {'C'}                                  % If the element is carbon, turn it into graphite, the standard form of carbon.
                cellarray{1,i}=[cellarray{1,i} '(gr) '];
            case {'S'}                                  % If the element is sulfur, turn it into the reference state for sulfur at that temp. 
                if T<368.3007
                    cellarray{1,i}=[cellarray{1,i} '(a)'];
                elseif T<388.3607
                    cellarray{1,i}=[cellarray{1,i} '(b)'];
                else
                    cellarray{1,i}=[cellarray{1,i} '(L)'];
                end
            otherwise
        end
    end
end