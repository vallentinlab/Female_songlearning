clc
close all
clear all

[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

%% load in workspaces:

load('MD_bird_stim_IDs_added23176.mat')

filelist=dir('*0ms_silence.mat')

for i=10%:length(filelist)

filename = char(strcat(pathname,'\',filelist(i,1).name));
load(filename)

number=size(traces,2);

%% demean traces:

for r = 1:number
averagemembranepotential = mean(traces{1,r}(1:end-1000,1),1);
demean_data(r,:) = traces{1,r}(:,1)-averagemembranepotential;;
traces_new{1,r}=demean_data(r,:).';
traces_demeaned{1,r}=traces_new{1,r};
end

 for k = 1:number

%% Remove spikes (Jagadeesh 1997) %%

spikes= find(traces_new{1,k}(:,1)>20);%-20);
for p = 1 : length(spikes)
        if spikes(p)<121 ||spikes(p)> size(traces_new{1,k}(:,1),1)-500
        continue
    end
traces_new{1,k}([spikes(p)-120:spikes(p)+120],1)=NaN;
end
end


for k = 1:number
  traces_new{1,k}(:,1)=fixgaps(traces_new{1,k}(:,1));
end

%% figure out which is BOS:

count=0;
count_fem=0;
for t=1:length(sound)
   
    figure(t)
    subplot(4,1,2:3)
    
    spectrogram((sound{1,t}),1024,1000,1024,40000,'yaxis');

colorbar off

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-60 5])
ylim([0 7])
set(gcf,'color',[1 1 1])

temp_idx(t)=bird.stim_ID{i}(t);

    if temp_idx(t)==1
        
         count=count+1;
         BOS_idx(count)=t;
    
         BOS_traces{1,count}=traces_new{1,t};
         
    elseif temp_idx(t)==0
        
        
        count_fem=count_fem+1;
        FEM_idx(count_fem)=t;
        
        FEM_traces{1,count_fem}=traces_new{1,t};


    end

close ALL
clear temp_idx

end

%% calculate the mean of BOS trace:

mat_BOS_traces=cell2mat(BOS_traces).';

if length(BOS_idx)>1

avg_BOS_trace=mean(mat_BOS_traces);

elseif length(BOS_idx)==1
    
    avg_BOS_trace=mat_BOS_traces;

end

song_envelope=smooth(abs(hilbert(sound{1,1})),50);


if i==5;
    
    FEM_idx(3)=FEM_idx(2);
    FEM_idx(2)=FEM_idx(1);
end

%% overlay traces for each data point:

for z=1:length(FEM_idx)

figure_handle=figure(z)

plot1=subplot(4,1,1);

 spectrogram((sound{1,BOS_idx(1)}),1024,1000,1024,40000,'yaxis');

colorbar off

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-60 5])
ylim([0 7])

plot2=subplot(4,1,2);

 spectrogram((sound{1,FEM_idx(z)}),1024,1000,1024,40000,'yaxis');

colorbar off

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-60 5])
ylim([0 7])

 
linkaxes([plot1, plot2], 'x' ); %link X-axes together in order to zoom in simultaniously in all plots

axes1=findobj(plot1, 'type', 'axes');
axes2=findobj(plot2, 'type', 'axes');
all_ha = findobj( figure_handle, 'type', 'axes');

%% figure out x and y values by clicking on the plot:
% first the program pauses in order for you to zoom in
% when you have zoomed in to onset, press any key (space)
% click left mouse click to set onset, then klick right mouse click, if you
% are done. Then zoom out and zoom into offset. Again press space, set
% coordinates, finish with right mouse click. Your x-coordinates will appear
% in variable 'timing'.

for u=1:2

pause

button=1;
while sum(button)<=1
    [x_value, y_value, button]=ginput(1);
    timing(u)=x_value;
end
end
 

if length(sound{1,1})/40000<1 % check if sound is longer than 1 sec, different spectrogram
   
call_on(z)=round(timing(1,1));
call_off(z)=round(timing(1,2));
else
call_on(z)=round(timing(1,1)*1000);
call_off(z)=round(timing(1,2)*1000);

end

FEM_timing{1,i}(z,1)=call_on(z);
FEM_timing{1,i}(z,2)=call_off(z);


clear timing

clear diff_BOS

clear diff call_on call_off diff_BOS

clear temp_diff_FEM diff_FEM_cell


end


 clear FEM_idx
 
clearvars -except diff_vec filelist pathname bird diff_BOS_vec...
    avg_diff_BOS_percell std_diff_BOS_percell avg_diff_FEM_percell...
    std_diff_FEM_percell Cell_FEM_max_coeff Cell_FEM_correlated_lag Cell_BOS_correlated_lag Cell_BOS_max_coeff...
    Cell_call_timing FEM_timing
end


%% plot correlated lags and xcorr coefficients:
figure1=figure(31);
for t=1:length(Cell_BOS_correlated_lag)
   
    color=rand(1,3);
    
    control_xcorr{1,t}=reshape(Cell_BOS_max_coeff{1,t},1,[]);
    treat_xcorr{1,t}=reshape(Cell_FEM_max_coeff{1,t},1,[]);
    
   ctrl_xcorr_mean(t)=mean(control_xcorr{1,t});
   treat_xcorr_mean(t)=mean(treat_xcorr{1,t});
    
    control{1,t}=reshape(Cell_BOS_correlated_lag{1,t},1,[]);
    treat{1,t}=reshape(Cell_FEM_correlated_lag{1,t},1,[]);
    
    ctrl_mean(t)=mean(control{1,t})/40; % make it in ms instead of units
    treat_mean(t)=mean(treat{1,t})/40; %same
    ctrl_std=std(control{1,t});
    treat_std=std(treat{1,t});
    
    ctrl_dots=repmat(1,1,length(control{1,t}));
    treat_dots=repmat(3,1,length(treat{1,t}));
    

    
    figure(40+t)
    
    plot(ctrl_dots,control{1,t},'o','Color', color, 'MarkerSize', 20)
    hold on
    plot(treat_dots,treat{1,t}, 'o','Color', color, 'MarkerSize', 20)
    hold on
    xlim([0 4])
    

end

x1=cell2mat(control);
x2=cell2mat(treat);

sem_xcorr_contr=std(ctrl_xcorr_mean)/sqrt(length(ctrl_xcorr_mean));

sem_contr=std(ctrl_mean)/sqrt(length(ctrl_mean));

stats=ranksum(ctrl_mean,treat_mean);

all_means=[ctrl_mean; treat_mean].';

all_xcorr_means=[ctrl_xcorr_mean; treat_xcorr_mean].';

figure1

subplot(1,2,2)

for f=1:length(all_means)

plot((1:2),all_means(f,1:2), '-o')
    hold on
plot([0 3], [mean(ctrl_mean)+1.95*sem_contr mean(ctrl_mean)+1.95*sem_contr], 'r-')
hold on
plot([0 3], [mean(ctrl_mean)-1.95*sem_contr mean(ctrl_mean)-1.95*sem_contr], 'r-')
xlim([0 3])
end
axis square
box off
ylabel('Lag time (ms)')
set(gca,'XTickLabel',{'','', 'Song snippet', '', 'Song+call', ''});
set(gcf,'color',[1 1 1])


subplot(1,2,1)

for f=1:length(all_xcorr_means)

plot((1:2),all_xcorr_means(f,1:2), '-o')
    hold on
plot([0 3], [mean(ctrl_xcorr_mean)+1.95*sem_xcorr_contr mean(ctrl_xcorr_mean)+1.95*sem_xcorr_contr], 'r-')
hold on
plot([0 3], [mean(ctrl_xcorr_mean)-1.95*sem_xcorr_contr mean(ctrl_xcorr_mean)-1.95*sem_xcorr_contr], 'r-')
xlim([0 3])
end
axis square
box off
ylabel('Max subthreshold correlation')
set(gca,'XTickLabel',{'','', 'Song snippet', '', 'Song+call', ''});
set(gcf,'color',[1 1 1])

stats=ranksum(ctrl_xcorr_mean,treat_xcorr_mean);
stats_lag=ranksum(ctrl_mean,treat_mean);
