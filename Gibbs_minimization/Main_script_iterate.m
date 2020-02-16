%%% This script calculates the available Gibbs energy for the specified
%%% multiphase system by repeatedly minimizing the total Gibbs energy of
%%% the system starting from several (random) initial conditions.

warning off all
vals=[];  % this array will store the observed and equilibrium abundances for each iteration
% in addition to derived properties like the Gibbs free energy difference between
% the intitial state and the equilibrium state, delta_G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations:
num=30; %~10 should give roughly the global minimum, ~100 will reliably provide an accurate global minimum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_vectors %this load the names, species, and abundances of species
n=n_init.*scale_factor; %multiply by large scale factor for more reliable minimization
% routine will divide abundances by scale_factor at the end to return to
% sensible units.

% iterate over differnt initial conditions
 for i=1:num
    i
    n=randperm(length(n));
    for j=1:length(n)
        n(j)=n(j)*500*rand()*10^(10*rand()-9); %this generates totally arbitrary random initial conditions for each iteration
    end
    tic
    vals=[vals; n Gibbs_energy_minimization(n')]; % run the Gibbs energy minimization routine and store the output in the vals array
end
n = n_init.*scale_factor

%find global minimum value
[M,I]=min(vals(:,length(n)*2+3));

% Diplsay Gibbs energy difference from all iterations
disp(['Array containing available Gibbs energies from each iteration:'])
array_of_Gibbs=vals(:,length(n)*2+3)/scale_factor

% re-load initial state vectors for plotting purposes
load_vectors;
n_true = n_init;
l=l_init';

 %%% Option to re-normalize mixing ratios to 1 (not recommended)
 % [row,ve] = find(l==1);
 % n_true(row)=n_true(row)/sum(n_true(row));


% change molalities to moles per mole of atmosphere for AQ species
 mass_ocean= om*1.3802e21; % mass of ocean in kg
 moles_atm=1.762035714285714e+20; % total number of moles in Earth's atmosphere
 [row,ve] = find(l==4); %% find indices for aqueous species only
 n_true(row)=n_true(row).*mass_ocean./moles_atm; % convert molaltities to moles per mole atmosphere

% Extract final abundances for the global minimum case:
min_val=length(n_true)+1;
max_val=min_val+length(n_true)-1;
final_n=vals(I,min_val:max_val)'/scale_factor; % final abundances global minimum case

G_dex=max_val+3; % index for Gibbs energy difference
disp(['Global minimum (J/mol):'])
deltaG_value=vals(I,G_dex)/scale_factor %display the largest (negative) Gibbs energy difference
if deltaG_value > 0 % display warnings if Gibbs energy change is positive
    disp(['Warning, postive Gibbs energy change - global minimum has not been obtained'])
    disp(['Try more iterations'])
end

%%% Plot outputs
%Plot_outputs(n_true',final_n)
toc
%%% save result
save_result
