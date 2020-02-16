%%% This function plots the initial and final states as histograms.
%%% Note that for systems other than the Archean and Proterozoic examples,
%%% this portion of this script that produces the second figure will 
%%% need to be modified to include the correct species.

function outs=Plot_outputs(n_true,final_n,names)

warning off all                     



% Create first figure with all species on same axis
% This should work for any input file
figure
cutoff=1e-10; %% Where to cutoff very small abundances
mycolor=[0.3 0.1 0.8;0.2 0.2 0.2;0.1 0.95 0.4];
colormap(mycolor) 
set(gca,'YScale','log','fontsize',10)
hold on;

str=names;
names
                                                     
bar([n_true,final_n],'basevalue',cutoff) 
set(gca,'XTicklabel',str,'Xtick',1:numel(str))
ylabel('Mixing ratio')
legend('observed','equilibrium')
title('Earth')
grid on


