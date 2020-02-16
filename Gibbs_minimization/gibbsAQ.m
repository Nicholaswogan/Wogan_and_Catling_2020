% This function returns the value of the Gibbs Free Energy of Formation 
% for aqueous species at the specified temperature and pressure.
% Throws errors if the designated species cannot be located in the
% database or if the species do not have an acceptable temperature range.

% This function uses equation 32 in Krissansen-Totton et al. (2016) to
% calculate Gibbs energies of formation. Parameters are taken from the
% sprons96 database.
function G=gibbsAQ(species1,T,P)
    global databaseC;                % Load the database
    species = regexptranslate('wildcard',species1); % Makes it so that any terms (like +) translate literally when searching.
    index = searchDataC(species);    % Get the index in the database of the species
    
    if index == -1                  % If the species isn't there, throw an exception
        missingspecies = species1
        throw(MException('Catling:missingSpecies',['Species ' species1 ' is missing from the database.']))
    end
    if index ==-2                   % If there is more than one match, throw an exception.
        notSpecificSpecies = species1
        throw(MException('Catling:nonSingularSpeciesName',['Species ' species1 ' must be more specific to find a match.']))
    end
    databaseC{index,1};
    coef = databaseC{index,2};       % Get the coefficient and temperature range data from the database.
      
    if isnumeric(coef)              % If the coefficients are in an acceptable form, calculate G.
        % See equation 32 in Krissansen-Totton et al. (2016) for a
        % definition of all the variables and coefficients below:
        Tr = 298.15;
        Pr = 1;
        Psi = 2600;
        Theta = 228;
        Y = -5.81*10^(-5);
        Gr = coef(1);
        Hr = coef(2);
        Sr = coef(3);
        a1 = coef(4)/10;
        a2 = coef(5)*100;
        a3 = coef(6);
        a4 = coef(7)*10000;
        c1 = coef(8);
        c2 = coef(9)*10000;
        w = coef(10)*100000;
        q = coef(11);
        diE = dielectric(T,P);
                            % These formulas are all parts of the formulas used in Walther's Essentials of GeoChemistry.
        G1 = Gr;
        G2 = -1*Sr*(T-Tr);
        G3 = -1*c1*(T*log(T/Tr)-T+Tr);
        G4 = a1*(P-Pr);
        
        h1 = log( (Psi+P) / (Psi + Pr));
        
        G5 = a2*h1;
        
        h2 = (1/(T-Theta))-(1/(Tr-Theta));
        h3 = (Theta-T)/Theta;
        h4 = log(Tr*(T-Theta)/(T*(Tr-Theta)));
        
        G6 = -1*c2*(h2*h3-T/(Theta*Theta)*h4);
        G7 = (1/(T-Theta))*(a3*(P-Pr)+a4*h1);
        G8 = w*Y*(T-Tr);
        
        G = 4.184*(G1+G2+G3+G4+G5+G6+G7+G8);
    end

    
    clear('index','data','coef','a','H','S') % Clear all temporary variables used.
end