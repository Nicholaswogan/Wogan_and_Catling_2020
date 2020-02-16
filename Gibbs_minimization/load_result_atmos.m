fileID = fopen('Volc_iter_7.8.19/Q25/Gibbs_result_pre_atmos_only.txt', 'r');
%fileID = fopen('result_1st_eco_middle.txt', 'r');

for i=1:3
    fgetl(fileID);
end



AAA =[];
%AAA = fscanf(fileID,'%f %f',[2 Inf]);
for i=1:9
    a = str2num(fgetl(fileID));
    AAA = [AAA; a];
end




%AAA = AAA';

n_input = AAA(:,1);
n_output = AAA(:,2);




fclose(fileID);

Plot_outputs_atmos(n_input,n_output,names)