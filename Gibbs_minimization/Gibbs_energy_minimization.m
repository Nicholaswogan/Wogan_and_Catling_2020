%%% The function Gibbs_energy_minimization takes a (random) initial state
%%% and calculates the multiphase equilibrium appyling atom and charge
%%% conservation constraints, which are derived from the true initial state
%%% (in load_vectors). It returns equilibrium abundances and the Gibbs energy difference
%%% between the (true) initial and equilibrium states. The only parameters that 
%%% should be modified in this file are temperature and pressure, and potentially
%%% some of the fmincon parmaters, if the optimization routine is having 
%%% trouble finding minimia.

function vals=Gibbs_energy_minimization(n)
  
    
    clear global R T P l v g V names;      % Clear any global variables.
    format('long')                         % Formats the data output to reveal more decimal places.
    global R T P l v g V names n_true a useful om scale_factor;            % Allows values set in this method to be passed to other methods.
    
    R = 8.3145;                            % Universal Gas Constant (J/K/mol)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = 298.00;            % Temperature of system in Kelvin
    P = 1.0;                 % Pressure of system in Bars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load input files (abundances, species names, charges etc.)
    load_vectors
    % names = cell array containing the names of each element in the simulation. 
    l=l_init'; % l = matrix containing the phase of each species.
    v=v_init;     % v = matrix containing the charge of each species.

     
    
    elements = makeAs(names);                      % Returns a cell array containing the '# of elements per mole of each species' matrix.
                                                   % Remove the element names at the top of the matrix to convert to a matrix.
    a = cell2mat(elements(2:end,:))';              % Coefficient Matrix a(ij) for element i, molecule j.
 

    % Loop to calcualte the gibbs free energy of each species.    
    for i=1:length(names) 
        if l(i) == 4                      % If the species are aqueous, use the aqueous database. Otherwise, use the gaseous database.
            g(i) = gibbsAQ(names{i},T,P); % GibbsAQ returns the aqueous Gibbs free energy values of formation
        else
            g(i) = gibbs(names{i},298.15)+gibbsBB(names{i},T)-gibbsBB(names{i},298.15); %% Gibbs energy of formation for gaseous species
            % This formulation is required to ensure the gaseous Gibbs energies of formation 
            % match the conventions for the aqueous species Gibbs energies of formation.
        end
    end
    

    n_true = n_init'; %% Use the true initial state to calculate atom and charge conservation conditions
    
    % option to force moles in atmosphere to sum to one (not recommended)
    %[row,ve] = find(l==1);
    %n_true(row)=n_true(row)/sum(n_true(row));
    
     % convert molalities to moles per mole atmosphere for AQ species
     mass_ocean=om*1.3802e21; % number of moles in ocean
     moles_atm=1.762035714285714e+20; % number of moles in atmosphere
     [row,ve] = find(l==4); % find indices for aqueous species
     n_true(row)=n_true(row).*mass_ocean./moles_atm; %convert to moles per mole atmosphere

     %ensure net charge is zero by adding Na(+)
     IndexC =strfind(names, 'Na(+)');
     Index = find(not(cellfun('isempty', IndexC)));
     n_true(Index)=n_true(Index)-v*n_true; %adjsut Na(+) to ensure zero net charge
    useful=n_true;
    
    n_true=n_true.*scale_factor; % Multiple abundances by scaling factor for more precise calculations at low abundances

    G1=totalGibbsInternal(n_true); % Calculate the total Gibbs energy of the true initial state.
   
     
     b=a*n_true; % atom conservation constraint that will be used in minimization routine.
     V=v*n_true; % charge conservation constraint (total charge should be close to zero) that will be used in minimization routine.

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADJUSTABLE VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %                These values can be modified to change the accuracy of the algorithm.
    tolfun = 1e-60;                                        % The required tolerance in the Gibbs free energy before exiting the minimization algorithm
    tolX = 1e-60;                                          % The required tolerance in n values before exiting the minimization algorithm
    tolMin = 1e-60;                                        % Smallest possible real number for the input parameters.
    tolCon = 1e-60;                                        % Maximum value of the constraint function  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Options, which can be set to modify the behavior of the fmincon minimization algorithm
    options = optimset('Display','on','GradConstr','off','GradObj','on','DerivativeCheck','off','FinDiffType','central','TolX',tolX,'TolFun',tolfun,'TolCon',tolCon,'FunValCheck','on','Algorithm','interior-point','AlwaysHonorConstraints','bounds','MaxFunEvals', 500000,'MaxIter', 10000); 
    % Run fminon minimization routine, starting from random initial state
    % (n) but using atom conservation constraints (a and b) and charge
    % conservation constraints (@mycon) from true initial state. n_true.
    [n,G2,exitflag,s]=fmincon(@totalGibbsInternal,n,[],[],a,b,repmat(tolMin,size(n)),Inf,@mycon,options);                                             % a  b
    % fmincon returns the equilibrium abundances (n), total Gibbs energy of
    % the final state (G2), and the exitflag for the minimization routine.
    
    exitflag %display exit flag for each iteration
    
    G2;              % The final total Gibbs free energy.
    dG = G2-G1;      % The change in Gibbs free energy of the system between initial and final state.
    
    % Output selected results.
    disp(['Mass Balance "Roughness" (should be zero):           ' num2str(sqrt((b-a*n)'*(b-a*n)))])
    disp(['Charge Balance (should be zero):                     ' num2str(V-v*n)])
    disp(['dG value:                                            ' num2str(dG)])
    disp(['G1 value:                                            ' num2str(G1)])  
    disp(['G2 value:                                            ' num2str(G2)])
    disp(['scaled dG value:                                     ' num2str(dG*8.3145*T/scale_factor)])

     %%% Fill vals array with selected outputs
     for zz=1:length(names)
         vals(zz)=n(zz);
     end
     vals(length(names)+1) = dG;
     vals(length(names)+2) = G2;
     if exitflag==0 %% Use fill values of 9999.99 if fmincon was unsuccessful
        vals(length(names)+3) = 9999.99*scale_factor;
     else
         vals(length(names)+3) = dG*R*T;
     end
     
      %%% Display initial and final abundances for every species
      for ii=1:length(names)
          [names(ii) n_true(ii)/scale_factor n(ii)/scale_factor]
      end

end

%%% The function totalGibbsInternal calculates the total Gibbs energy of
%%% the system given abundances, n, and the temperature and pressure of the
%%% system.
function [fun,grad] = totalGibbsInternal(n)
    global R T g P l names v n_true scale_factor om;         % Loads the variables from the above function.
    
    TotalMolesGases = sum(onlyPhase(1,n));
    TotalMolesLiquids = sum(onlyPhase(0,n))+sum(onlyPhase(2,n))+sum(onlyPhase(4,n));
    TotalMolesSolids = sum(onlyPhase(3,n));
         
    % The natural log of the fugacity coefficient for each of the gas phase species in the system
    lnPhi = fugCoef(T,P,onlyPhase(1,names),onlyPhase(1,n'));
    lnPhiAQ = fugCoefAQ(T,P,onlyPhase(4,names),onlyPhase(4,n'),onlyPhase(4,v)); 
    %%% fugCoefAQ outputs are superceded by Pitzer_activity_diseq outputs (see below),
    %%% but the function is retained incase the user wishes to explore
    %%% regimes where the Truesdell-Jones equation is more appropriate than
    %%% Pitzer equations.
    
    %Calculate water activity and aqueous species activity coefficients
    %using the Pitzer equations
    [lnPhiWATER,act_an,act_ct,a_index,c_index]=Pitzer_activity_diseq(names,n,v,l);
    diff_for_index_conv=length(l)-length(lnPhiAQ);
    lnPhiAQ(a_index-diff_for_index_conv)=log(act_an); %% Fill in activity coefficients for anions
    lnPhiAQ(c_index-diff_for_index_conv)=log(act_ct); %% Fill in activity coefficients for cations
    
    % The contribution to the total gibbs free energy from all the species of each phase of the system.
    % This is essentially equation 10 in Krissansen-Totton et al. (2018):
    funct{1} = onlyPhase(0,n)'*(onlyPhase(0,g)/R/T+lnPhiWATER); %% 
    funct{2} = onlyPhase(1,n)'*(onlyPhase(1,g)/R/T+log(P)+lnPhi+log(onlyPhase(1,n)/TotalMolesGases));  %% 
    funct{3} = onlyPhase(2,n)'*(onlyPhase(2,g)/R/T+log(P)+log(onlyPhase(2,n)/TotalMolesLiquids));  %
    funct{4} = onlyPhase(3,n)'*(onlyPhase(3,g)/R/T+log(P)+log(onlyPhase(3,n)/TotalMolesSolids));  %
    funct{5} = onlyPhase(4,n)'*(onlyPhase(4,g)/R/T+lnPhiAQ+log(55.508435)+log(onlyPhase(4,n)/TotalMolesLiquids)-log(onlyPhase(0,n)/TotalMolesLiquids));  %%
   
    % Here, the contribution by each phase are totaled together into a total Gibbs free energy of the system.
    fun = 0;
    for i = 1:5                     % For each phase of matter.
       if ~isempty(funct{i})        % If the phase had any contribution, add it to the total gibbs free energy.
           fun = fun + funct{i};
       end
    end
   
    
    % Calculate analytic gradient of the system for each species in the system:
    % This follows equation 36 in Krissansen-Totton et al. (2016).
    grad = zeros(1,length(n));
    for i = 1:length(n)
        switch l(i)               % For each species, identify it's phase of matter and use the appropriate formula.
            case 0
                grad(i) = g(i)/R/T+log(n(i)/TotalMolesLiquids)-TotalMolesLiquids/n(i) - n(i)/TotalMolesLiquids + 2;
            case 1
                grad(i) = g(i)/R/T+log(P)+ lnPhi(1)+log(n(i)/TotalMolesGases); 
                lnPhi = lnPhi(2:end);  %  This removes the first fugacity coefficient from the list of coefficients. That way, each
                                       %       coefficient is only used once.
            case 2
                grad(i) = g(i)/R/T+log(P)+log(n(i)/TotalMolesLiquids);
            case 3
                grad(i) = g(i)/R/T+log(P)+log(n(i)/TotalMolesSolids);
            case 4
                grad(i) = g(i)/R/T + log(55.508435) + lnPhiAQ(1)+ log(n(i)/TotalMolesLiquids) - log(n(1)/TotalMolesLiquids) - n(1)/TotalMolesLiquids + 1; 
                lnPhiAQ = lnPhiAQ(2:end);
               
        end
    end 
    
end


% This function returns a form of some original matrix (old) with only species in it of the specified phase (phaseNumber)
function new = onlyPhase(phaseNumber, old)           %Takes a given matrix and removes any values from that matrix that are not associated with.
    global l                                         % the specificed phase.
    new = [];                                        % The new version of the matrix with only the given phase.
    for temp = 1:length(old)                         % For each species,
        if l(temp)==phaseNumber                      %  if the species phase is the correct phase
            new = [new;old(temp)];                   %  add it to the new matrix (thus, the new matrix will only contain species of the
        end                                          %  correct phase.
    end
end

% This function is used by fmincon to represent the charge balanace.
function [c,ceq,cd, ceqd] = mycon(x)
    global v V; % v = charges for each molecule; V = total charge of the system (should be 0 with sensible inputs)
    c = [];     % Compute nonlinear inequalities at x. This is blank since we have no inequalities
    ceq = (v*x-V);   % Compute nonlinear equalities at x (the moles values).
    cd = [];
    ceqd = v';
end







