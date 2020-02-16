fileID = fopen('Volc_iter_7.8.19_water_vapor/Q1/Gibbs_result_pre.txt', 'r');
%fileID = fopen('result_1st_eco_middle.txt', 'r');

for i=1:3
    fgetl(fileID);
end



AAA =[];
%AAA = fscanf(fileID,'%f %f',[2 Inf]);
for i=1:26
    a = str2num(fgetl(fileID));
    AAA = [AAA; a];
end




%AAA = AAA';

n_input = AAA(:,1);
n_output = AAA(:,2);




fclose(fileID);
close all
Plot_outputs(n_input,n_output)