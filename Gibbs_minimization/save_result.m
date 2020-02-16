fileID = fopen('Gibbs_min_result.txt', 'w');

AAA = [n_true' final_n];

% Save Gibbs energy value
fprintf(fileID,'Gibbs Energy (J/mol):\n');
fprintf(fileID,'%6.4e',deltaG_value);

% Save initial and final values
fprintf(fileID,'\ninitial_n final_n\n');
fprintf(fileID,'%6.4e %6.4e\n', AAA');

fclose(fileID);