%%
% Calculate subthreshold DELTAS 
%

clear all
close all
clc
%% Load data/ adress folder %% 
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

load('Subthreshold_juvenline_HVC_deltas_silence_E.mat')

contr_mean_trace=mean_call_trace;

contr_list=filelist2;

clearvars -except contr_mean_trace contr_list

load('Subthreshold_juvenline_HVC_deltas_E.mat')

treat_mean_trace=mean_call_trace;

treat_list=filelist2;

clearvars -except contr_mean_trace contr_list treat_mean_trace treat_list

load('Significant_cell_values_subthreshold_E.mat')



for g=1:length(contr_mean_trace)
    
    temp_avg_contr=mean(contr_mean_trace{g,1}); % demean the trace first
    
    temp_contr_trace=contr_mean_trace{g,1}-temp_avg_contr;
    
    temp_treat_trace=treat_mean_trace{g,1}-temp_avg_contr; % demean it using the control
    
    [max_contr, id_contr]=max(abs(temp_contr_trace));
    
    sign_contr=sign(temp_contr_trace(id_contr));
    
    [max_treat, id_treat]=max(abs(temp_treat_trace));
    
    sign_treat=sign(temp_treat_trace(id_treat));
 

    Delta(g)=(max_treat*sign_treat)-(max_contr*sign_contr);
    
    clear temp_countr_max temp_treat_max idx_contr idx_treat temp_treat_trace temp_contr_trace...
        temp_avg_contr temp_avg_treat max_contr max_treat id_treat id_contr sign_treat sign_contr...
    
end

figure

a=repmat(1,length(cell_idx),1);

scatter(a,Delta(cell_idx),'mo','Jitter', 'on')
hold on
plot([0 2], [0 0], '--k', 'LineWidth',2)
axis square
box off
xlim([0 2])
ylim([-10 10])