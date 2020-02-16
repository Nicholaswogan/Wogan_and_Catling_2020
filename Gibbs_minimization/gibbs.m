% This function returns the value of the Gibbs Free Energy of Formation at
% 1 ATM for the given gaseous species at the specified temperature.
% Throws errors if the designated species cannot be located in the
% database or the species does not have an acceptable temperature range.


% This shell function makes it so that the function can accept either a cell array of multiple species, or a single species.
function G=gibbs(species,T)
    a = makeAs(species);    % Calculates the coefficient matrix 'a' (ie: CH4 is C 1 H 4)

    if (iscellstr(species)) % If its a list, then perform the function on them all.
        for i=1:length(species)  % For each species
            G(i)=gibbsHelper(species(i),T); % Calculate the apparent gibbs free energy of that species
        end
        G=G'; % Flip the result so that the matrix algebra works later.
    else
        G=gibbsHelper(species,T);   %If it is only one species, just calculate the gibbs free energy of the specific species.
    end

    

    %%%% For converting Berman-Brown convention to Standard Gibbs free energies of formation
    elementSpecies = convert2standard(a(1,:),T);   % Obtain a list of all of the reference state species for each element.
    speciala = makeAs(elementSpecies);            % Create a cell array containing the number of each reference state species is in each species
    elementg=[];                                   % The gibbs free energy contribution by the elements
    for i=1:length(elementSpecies)                 % For each elemental species.
        ig = gibbsHelper(elementSpecies{i},T);     %    calculate the gibbs free energy contribution by that species.
        ig = ig/speciala{1+i,i};                  %    adjust that value by the number of moles each element has in it (H2 has 2, for instance)
        elementg = [elementg; ig];               %    add that new value to a matrix.
    end
    
    G = G - cell2mat(a(2:end,:))*elementg;         % With all of those matricies, multiply those values by the number of each element the
                                                % species have in them. (So, since CH4 has 4 H's, subtract the gibbs free energy of
                                                  % formation value by 4 * (g_H2)/2
end

function G=gibbsHelper(species,T)
    global database;                % Load the database
    species = regexptranslate('wildcard',species); % Makes it so that any terms (like +) translate literally when searching.
    index = searchData(species);    % Get the index in the database of the species
    
    if index == -1                  % If the species isn't there, throw an exception
        missingspecies = species
        throw(MException('Catling:missingSpecies',['Species ' species ' is missing from the database.']))
    end
    if index ==-2                   % If there is more than one match, throw an exception.
        notSpecificSpecies = species
        throw(MException('Catling:nonSingularSpeciesName',['Species ' species(:) ' must be more specific to find a match.']))
    end
    
    data = database{index,2};       % Get the coefficient and temperature range data from the database.
    coef = [];                      % Declare the coefficient matrix to be used for calculating G.
    
    for a = 1:min(size(data))       % Find the correct temperature range and store that set of coefficients into the variable 'coef'.
        if (data(1,a)<=T)&&(data(2,a)>=T)
            coef = data(3:end,a);
            break
        end
    end
    
    %if length(coef)<2               % If none of the temperature ranges fit the designated temperature, throw an exception.
    %    throw(MException('Catling:outsideTemperatureBounds',['Species ' species ' does not have coefficients for the specified temperature range.']))
    %end
    
    if length(coef)<2              % If none of the temperature ranges fit the designated temperature
        %disp(['WARNING: Species ' species ' does not have coefficients for the specified temperature range.'])
        if T < data(1,1)         % If the temperature is less than the lowest temp
            coef = data(3:end,1);    % Assign coef to be the lowest temperature range possible (which will hopefully be close).
        else
            coef = data(3:end,end);  % Otherwise, use the highest range (since that should be closest to the higher temperatire).
        end
    end
    
    coef;
    
    if isnumeric(coef)              % If the coefficients are in an acceptable form, calculate G.
        H = -1*coef(1)*T^(-2) + coef(2)*log(T)/T + coef(3) + coef(4)*T/2 + coef(5)*T^2/3 + coef(6)*T^3/4 + coef(7)*T^4/5 + coef(9)/T;
        S = -1*coef(1)*T^-2/2 - coef(2)*T^-1 + coef(3)*log(T) + coef(4)*T + coef(5)*T^2/2 + coef(6)*T^3/3 + coef(7)*T^4/4 + coef(10);
        G = 8.3145*T * (H - S);
    end
    %newH=H*8.3145*T;
    %newH
    clear('index','data','coef','a','H','S') % Clear all temporary variables used.
end