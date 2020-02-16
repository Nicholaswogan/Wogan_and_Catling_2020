%%% This script loads the input text file containin species names,
%%% abundances, charges, and species type, and converts these to vectors
%%% and cells that can be used by the main code.

global input_file

% Choose input file (uncomment the input you wish to use, or create your own)
%input_file = 'inputs_Archean_max.txt';
%input_file = 'inputs_Archean_min.txt';
%input_file = 'inputs_Proterozoic_max.txt';
%input_file = 'inputs_Proterozoic_min.txt';

input_file = 'inputs_Volc_iter.txt';
%input_file = 'inputs_Volc_iter_atmos_only.txt';

All_inputs=importdata(input_file,';'); %load from text file
v_input=All_inputs(1); % first line is charge vector
l_input=All_inputs(2); % second line is phase vector
n_input=All_inputs(4); % fourth line is abundance vector
scalar_input=All_inputs(5); % fifth line contains various scalars
name_input=horzcat('names = {',cell2mat(cellstr(All_inputs(3))),'};'); % 3rd line contains names of species
eval(name_input); % load species names
v_init=str2num(cell2mat(v_input)); %convert charge cell to matlab vector
l_init=str2num(cell2mat(l_input)); %convert phase cell to matlab vector
n_init=str2num(cell2mat(n_input)); %convert abundance cell to matlab vector
scalars_in=str2num(cell2mat(scalar_input)); %convert scalar cell to matlab vector
om=scalars_in(1); % defines number of ocean masses
sms=scalars_in(2); % salinity scalar (no longer used)
scale_factor=scalars_in(3); % scale factor
n_init(1)=om*n_init(1); % multiply liquid water abundance by number of ocean masses
