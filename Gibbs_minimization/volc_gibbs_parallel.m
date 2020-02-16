%%% This script runnns the Gibbs solver in parallel for all of the
%%% volcanism iteration cases.


tic
QQ = 1:25;

parfor Q=QQ
    Main_script_iterate_func('inputs_Volc_iter.txt','test.txt')
    
end
toc