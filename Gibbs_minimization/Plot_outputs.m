%%% This function plots the initial and final states as histograms.
%%% Note that for systems other than the Archean and Proterozoic examples,
%%% this portion of this script that produces the second figure will 
%%% need to be modified to include the correct species.

function outs=Plot_outputs(n_true,final_n)
global input_file names
warning off all                     

%Determine if plotting Archean or Proterozoic results:
TP=isempty(strfind(input_file,'Proterozoic'));
TA=isempty(strfind(input_file,'Archean'));
TA = 0;

Eon='None'; %% This will contain Proterozoic or Archean string
if TP == 0
    Eon = 'Proterozoic';
elseif TA == 0
    Eon = 'Archean';
end

% Create first figure with all species on same axis
% This should work for any input file
figure
cutoff=1e-10; %% Where to cutoff very small abundances
mycolor=[0.3 0.1 0.8;0.2 0.2 0.2;0.1 0.95 0.4];
colormap(mycolor) 
set(gca,'YScale','log','fontsize',10)
hold on;

str=names;
                                                     
bar([n_true,final_n],'basevalue',cutoff) 
set(gca,'XTicklabel',str,'Xtick',1:numel(str))
ylabel('Mixing ratio')
legend('observed','equilibrium')
title('Earth')
grid on

% Create second figure with separate subplots for aqueous and gaseous 
% species. This may need to be edited for inputs different to the default
% Archean and Proterozoic examples.
h = figure('units','centimeters','position',[5 5 43 22])
subplot(2,1,1)
cutoff=1e-10; % Where to cutoff very small abundances: 7 for Archean max, 8 archean min
for i = 1 : length(n_true)
    if n_true(i)<cutoff
        n_true(i)=cutoff;
    end
    if final_n(i)<cutoff
        final_n(i)=cutoff;
    end
end
colormap(mycolor) 
set(gca,'YScale','log','fontsize',22)
hold on; 
if strcmp(Eon,'Archean')
    str={'O_2','N_2', 'H_{2}O','CO_2','NH_3','CH_4','CO','H_2','H_{2}S'}; %Archean
    bar([n_true(2:10),final_n(2:10)],'basevalue',cutoff) % Archean
elseif strcmp(Eon,'Proterozoic')
    str={'O_2','N_2','H_{2}O','CO_2','NH_3','CH_4','H_2','N_{2}O'}; %Proterozoic
    bar([n_true(2:9),final_n(2:9)],'basevalue',cutoff) % Proterozoic 
else
    print('Neither Proterozoic nor Archean - errors may be present in plotting')  
end

ylabel('Mixing ratio')
legend('initial','equilibrium')
title('(a) gaseous species','Fontsize',22,'FontWeight','normal')
ylim([cutoff 1])
set(gca, 'YTick', [1e-10 1e-8 1e-6 1e-4 1e-2 1],'XTicklabel',str,'Xtick',1:numel(str))
set(gca,'XMinorTick','off','YMinorTick','off')
%[hx,hy]=format_ticks_v2(gca,str);
grid on
grid minor
grid minor

if strcmp(Eon,'Archean')
    init_aq=[n_true(1)',n_true(11:26)']'; %Archean
    final_aq=[final_n(1)',final_n(11:26)']'; %Archean
    str ={'H_{2}O_{(L)}','H^+','Na^+','Cl^-','HCO_{3}^-','CO_{2}','CO_{3}^{2-}','OH^-','NH_{3}','NH_{4}^+','N_{2}','O_{2}','CH_{4}','CO','H_{2}S','SO_{4}^{2-}','H_{2}'}; %Archean
elseif strcmp(Eon,'Proterozoic')
    init_aq=[n_true(1)',n_true(10:25)']'; %Proterozoic
    final_aq=[final_n(1)',final_n(10:25)']'; %Proterozoic
    str ={'H_{2}O_{(L)}','NO_{3}^-','H^+','Na^+','Cl^-','SO_{4}^{2-}','HCO_{3}^-','CO_{2}','CO_{3}^{2-}','OH^-','NH_{3}','NH_{4}^+','N_{2}','O_{2}','CH_{4}','H_{2}S','SO_2'}; %Proterozoic
else
    print('Neither Proterozoic nor Archean - errors may be present in plotting')  
end


subplot(2,1,2)
cutoff=1e-8; % Where to cutoff very small abundances: 7 for Archean max, 8 for Archean  min 
ylim([cutoff 500])
colormap(mycolor) 
set(gca,'YScale','log','fontsize',22)
hold on;
bar([init_aq,final_aq],'basevalue',cutoff) 
set(gca,'XTicklabel',str,'Xtick',1:numel(str))
[hx,hy]=format_ticks_v2(gca,str);
y_string=['Moles per mole',char(10),'of atmosphere'];
ylabel(y_string)
legend('initial','equilibrium')
title('(b) aqueous species','Fontsize',22,'FontWeight','normal')
set(gca, 'YTick', [1e-8 1e-6 1e-4 1e-2 1 100]) % 8 for min Archean 6 for max Archean
set(gca,'XMinorTick','off','YMinorTick','off')
grid on
grid minor
grid minor


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MySavedFile','-dpdf')

end


