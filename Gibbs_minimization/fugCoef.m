% This function returns the value of the fugacity coefficient at
% the given temperature and pressure for the set of gaseous species
% with the specified mole amounts. This uses equations 27 and 28 in
% Krissansen-Totton et al. (2016).
% Throws errors if the designated species cannot be located in the
% database or if the computer has trouble solving the cubic equation
% for the compressability factor.

function lnPhi=fugCoef(temperature,pressure,names,n)
    global databaseB;                 % Load the database
    
    % Pre-alocate Space for the following variables
    index = zeros(length(names),1);   % The location of each species in the database.
    data = zeros(length(names),3);    % The critical temperature, pressure and acentric factor data for each species
    ai = zeros(length(names),1);      % The coefficient values of 'a' for each species.
    bi = zeros(length(names),1);      % The coefficient values of 'b' for each species
    alphai = zeros(length(names),1);  % The coefficient values of 'alpha' for each species.
    
    R = 8.314472;                     % Ideal Gas Constant
   
    k = zeros(length(names));
    
    
    q = size(n);
    if q(1)>q(2)
        n=n';
    end
    
    for q = 1:length(names)                  % For each species passed.
        index(q) = searchDataB(names(q));    % Get the index in the database of the species
    
        if index(q) == -1                    % If the species isn't there, throw an exception.
            throw(MException('Catling:missingSpecies',['Species ' names{q} ' is missing from the database.']))
        end
        %disp(['Using species: ' names(q)])   % Display which species are being used.
        
        for p = 1:3                          % For each of the three types of data . . .
            temp = textscan(databaseB{index(q),p+1},'%f');
            data(q,p) = temp{1};             % Read and store the data for species q into slot p.
        end
                                            % Place each type of data into
                                            % an easier to read variable.
                                            
        tCrit = data(q,1);                   % The critical temperature in Kelvin 
        pCrit = data(q,2)*10;                % The critical pressure in bars (the data in the database has units of MPa, so *10 to convert
        ace   = data(q,3);                   % The acentric value for the species
        
        
        ai(q) = .42747*(R^2)*(tCrit^2)/pCrit;       % Calculate the coefficient values for the species
        bi(q) = .08664*R*tCrit/pCrit;
        alphai(q) = (1+(0.48508+1.55171*ace-0.15613*(ace^2))*(1-sqrt(temperature/tCrit)))^2;
        
    end
    
    aaTotal = 0;                            % Initial the variable representing coefficient 'a alpha, total'.
                                            % Take note of the difference between this and 'a alpha, w/ indicies'. 
                                            % They are different.
                                            
    for q = 1:length(ai)                    % Over the length of a and alpha . . .twice
        for p = 1:length(ai)
            aai(q,p) = (1-k(q,p))*sqrt(ai(q)*alphai(q)*ai(p)*alphai(p));    % Calculate the coefficient representing 'a alpha, w/ indicies'
            aaTotal = aaTotal + n(q)*n(p)/sum(n)^2*aai(q,p);     % Total the values to find the value of 'a alpha, total'
        end
    end


    
    bTotal = (n/sum(n))*bi;                               % Calculate the value of bTotal, the coefficient representing 'b, w/o indicies'
    A = aaTotal*pressure/((R*temperature)^2);          % Calculate the value of the coefficient representing 'A'
    BTotal = bTotal*pressure/(R*temperature);         % Calculate the value of the coefficient representing 'B, w/o indicies'
    Bi = bi*pressure/(R*temperature);                  % Calculate the value of the coefficient representing 'B, w/ indicies'
    cubic = [1 -1 (A-BTotal-BTotal^2) -A*BTotal];         % Coefficient terms to solve for the compressability factor, z.
    zi = roots(cubic);                                    % Use the customized method 'solveCubic' to calculate the value of z, 
                                                          % the compressability factor.
    
                                                          % roots returns all three roots, so the real root still needs to be solved
                                                          % for.

    z = Inf;                                              % Initialize the compressability factor to infinity to help solve for the real value
    
    for q = 1:3                                                          % For each root of z . . .
        if isreal(zi(q)) && zi(q)>0 && abs(zi(q)-1)<abs(z-1)             % Find the real root that is positive (z must be positive by definition) and
            z = zi(q);                                                   % is the closest value to 1.
        end                                                      
    end
    if z == Inf                                           % If no values of z fit that descritpion, something is wrong.
        throw(MException('Catling:badIncompressability','The value of z is either non-real or less than zero.'))
    end      
    % All coefficients have been found. Solve for the ln(Phi)
    
    PhiA = Bi/BTotal*(z-1);                               % These are all just smaller parts of the whole thing
    PhiB = log(z-BTotal);
    PhiC = A/BTotal*(Bi/BTotal-2/aaTotal*(aai*(n/sum(n))'));
    PhiD = log(1+BTotal/z);
    lnPhi = PhiA-PhiB+PhiC*PhiD;
    

    Phi = exp(lnPhi);

end