
clear all
close all
clc
%% Load data/ adress folder %%
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

%% choose .mat file that includes aligned sound and recording traces



% for silence right before call treat is 'sil_on' (if call is before song)

filelist2 = dir('*_E.mat');
treat=('sil_on'); % set to 'no' when analyzing call only stimulus and to sil_no to measure silence before

zoom=input('Would you like to look at plots per cell? yes(1)/no(2)');
nozerotrials=input('Would you like to include zero spike trials in calculating precision? yes(1)/no(2)')

number1 = length(filelist2);

%% load in the onset/offset library:

load('Solving_percentage_problem_onset_offset_library.mat'); % for head-fixed juveniles


fields=fieldnames(library);


%% variables:

fs = 40000; % Sampling rate

pre_silence = 1000; 
post_silence = 1000;

if strcmp(treat,'no')==1 % if looking at only one call, adjust silent period
   clear pre_silence post_silence
   
    pre_silence=500;
    post_silence=500;
    
elseif strcmp(treat, 'sil_no')==1
       clear pre_silence post_silence
   
    pre_silence=50;
    post_silence=50;

end

pretime = pre_silence*40;
posttime=post_silence*40;
i= 1;

bird=[];
bird_ID = [];
treatment = [];
stdv_fr_call=[];
stdv_fr_playback=[];
dph=[];
cell_ID=[];



for z = 1:7
    
filename = char(strcat(pathname11,'\',filelist2(z,1).name));
load(filename)

%% for juv data:

bird = (filename(end-27:end-24));
bird_ID = [bird_ID; bird];
cell = (filename(end-10:end-6));
cell_ID = [cell_ID; cell];
old_traces=traces;
thresh=15;
threshold=15;


shift = 20; %40
thr = 0.5;
window = 20; % 20 ms window for shuffling
pairs_new=[];
dist_sharp=20;
dist_ex = 0.5;
nA = 0.01;
all_spikes_psth=[];
psthdata=[];
ISI=[];
ISI_silence=[];
number_of_spikes=[];

stim = (length(traces{1,1})-pretime-post_silence*40)/40;
stim2 = (length(traces{1,1})-pretime-posttime); % stim length in digits
stim_length = fix(stim); % in ms
stim_length_sec =  stim_length/1000;

number = length(traces);


%% check if celll activity is above 1 Hz
for k=1:length(traces)

mVS2=traces{1,k}(:,1)-smooth(traces{1,k}(:,1),fs/80);%high pass filter it. here it is at 2ms 



[pksp,time_spike_playback] = findpeaks(mVS2,'MinPeakHeight',thresh,'MinPeakDistance', 80); % find spikes in the whole playback
number_of_spikes_playback(k) = length(findpeaks(mVS2,'MinPeakHeight',thresh,'MinPeakDistance', 80));

number_of_spikes_playback2(k) = length(time_spike_playback);

fr_playback(k) = abs(number_of_spikes_playback(k)/stim_length_sec);

end

fr_cell=mean(fr_playback);
fr_cell_all(z)=fr_cell; % calculate average firing rate of the cell across all playbacks


clear mVS2 thresh pksp time_spike_playback number_of_spikes_playback fr_playback fr_cell


if fr_cell_all(z)>=1 %% check if this cell has at least 1 Hz firing rate
    

%% detect where call is:

for h=1:size(fields,1) %because 5 birds

    if bird(end-2:end)==fields{h}(end-2:end)
        
        if strcmp(treat,'sil_on')==1
    
    call_onset_ms=library.(fields{h}).onset.(treat)-100 % take a 50 ms gap between silence before offset and call to avoid fake precision
    offset_call_ms=library.(fields{h}).offset.(treat)-100
    
        elseif strcmp(treat,'no')==1
            
            call_onset_ms=501;
            offset_call_ms=(length(traces{1,1})/40)-500;
            
        elseif strcmp(treat,'sil_no')==1
            
            call_onset_ms=250;
            offset_call_ms=349; 
            
        else
    
    call_onset_ms=library.(fields{h}).onset.(treat) 
    offset_call_ms=library.(fields{h}).offset.(treat)
        end
    end
end

onset_call_digit = call_onset_ms*40;
offset_call_digit = offset_call_ms*40+2000;%added 50 ms

call_length = fix(offset_call_digit-onset_call_digit); % digits
call_length_s = call_length/40000;
call_length_ms = call_length_s*1000;


call_matrix = zeros(number, call_length*40);
silence_matrix=zeros(number,pretime);

indReport=100:100:1000;

%% silence and playback traces


for k = 1:number
   call_traces{1,k}=traces{1,k}(onset_call_digit:offset_call_digit,1);
   call_traces{1,k}(:,1)=call_traces{1,k}(:,1)-smooth(call_traces{1,k}(:,1),fs/80);
end

%% De-mean data by substracting the mean value of each trace and plot raw traces%%


matrix_call = cell2mat(call_traces);
for k = 1:number
    averagemembranepotential_call = mean(cell2mat(call_traces));
    demean_data_call(:,k) = matrix_call(:,k) - averagemembranepotential_call(k);
end

matrix_call= demean_data_call'; 

if zoom==1

subplot(3,2,1)
for u=1:number
    plot(1/40:1/40:(length(matrix_call))/40,u*dist_sharp+[matrix_call(u,:)],'k','LineWidth',1)
    hold on
    xlim([0 length(matrix_call)/40])
    xlabel('Time (ms)')
    ylabel('de-meaned Membrane Potential (ms)')
    axis square
end

box off
end

win = 25;
edges = [0:win : length(matrix_call)];
psth = zeros(length(edges),1)';
spikes_check={};
subplot(3,2,3)

%% change here to look at playback or call
for i = 1:number
    spikes = find_spikes(matrix_call(i,:),threshold);
    spikes_check{1,i}=spikes;
    
%    
    dot = zeros(1,length(matrix_call)); % create a zeros vector at the lenght of call
    dot(spikes) = 1; % put a 1 for each detected spike at the zeros vector
    
    psthdata = [psthdata;dot];
     
    if zoom==1
    
        for j=1:length(spikes)
        line([spikes(j)/40 spikes(j)/40],[i-1 i],'Color','k','LineWidth',1)
        xlim([0 length(matrix_call)/40])
        xlabel('Time (ms)')
        ylabel('Trials')
        axis square
        
    end
    
    end
    %end
   
        refractory = diff(spikes/40); % calculates latencies between spikes [X(2)-X(1) X(3)-X(2)]
    %ref_silence = diff(spikes_silence/40);
    
    if refractory >1
        ISI =[ISI refractory];
    end
%     if ref_silence >1
%         ISI_silence =[ISI_silence ref_silence];
%     end
    number_of_spikes = [number_of_spikes length(spikes)];
end

%% use this section to further only look at trials with spikes (during call)

if nozerotrials==2

no_spike_trials=find(number_of_spikes==0); % find trials with no spikes

psthdata(no_spike_trials,:)=[];

number_of_spikes(no_spike_trials)=[]; % get rid of trials with not spikes

number=length(number_of_spikes);

end
%%

test_spike=find(number_of_spikes~=0); %test if spikes only occur in one row:

    if length(test_spike)>1



%     %% Peri- Stimulus Time Histogram
    psth = psth + histc(spikes/40,edges); % counts number of values in spikes that fall between in elements in the edges vector

    binning = 40;
[idx_spikes_row idx_spikes_col]= find(psthdata==1);
sorted= sortrows([idx_spikes_row idx_spikes_col]);
idx_spikes_row =sorted(:,1);
idx_spikes_col = sorted(:,2);



if length(idx_spikes_col) > 1 %%%% trying to analyse if there are any spikes in data


for i=1:number

    idx=find(idx_spikes_row == i);
    bin_spikes(i,:) = histc(idx_spikes_col(idx),[0:binning:size(psthdata,2)]); % sort spikes in bins
end

latencies=[];


for i = 1:number
    spiketimes_trigger = find(bin_spikes(i,:)==1); %finds spikes (real data) across call length of traces
    if length(spiketimes_trigger)==0
        continue
    end
   

        [x,spiketimes_test] = find(bin_spikes([1:i-1 i+1:end],:)==1);% find the other spikes in other trials
        
        for m=1:length(spiketimes_trigger)
            latency =spiketimes_test-spiketimes_trigger(m);
            latencies = [latencies; latency];
        end
        

end

if zoom==1

subplot(3,2,5)
    intervall = [-size(bin_spikes,2):10:size(bin_spikes,2)];
histogram = histc(latencies,intervall);
bar(intervall,histogram,'k');
xlabel('Latency (ms)')
ylabel('Number of spikes')
xlim([-size(bin_spikes,2) size(bin_spikes,2)])
axis square
box off

end

STA_number = sum(histc(sort(latencies), -window:1:window))/length(latencies);
shuffle_time = 1000;
for u=1:shuffle_time
  
    choose_nr_spikes = randsample(number_of_spikes,number,true); % choose as many spikes  as the traces from real data randomly with replacement (if more traces than spikes, inlcude repeated values in sample)
    shuffled_data = zeros(number, size(bin_spikes,2));

    for j=1:number
        if choose_nr_spikes(j)==1
            shuffled_data(j,ceil(size(bin_spikes,2)*rand(1)))=1;
        elseif choose_nr_spikes(j)==0
            continue
        else
 position = floor(size(bin_spikes,2)*rand(1,choose_nr_spikes(j)))+1; % assign the chosen nr of spikes new positions
 shuffled_data(j,position)=1; % put the new random spike position in tall-trial dot-raster matrix (spike=1)
        end
    end
    
    
    %% STA % spike triggered average

latencies_shuffled=[];


for i = 1:number
    spiketimes_trigger_shuffled = find(shuffled_data(i,:)==1);
    if length(spiketimes_trigger_shuffled)==0
        continue
    end
   

        [x,spiketimes_test_shuffled] = find(shuffled_data([1:i-1 i+1:end],:)==1);
        
        for m=1:length(spiketimes_trigger_shuffled)
            latency_shuffled =spiketimes_test_shuffled-spiketimes_trigger_shuffled(m);
            latencies_shuffled = [latencies_shuffled; latency_shuffled];
        end
        
    STA_number_shuffled(u) = sum(histc(sort(latencies_shuffled), -window:1:window))/length(latencies_shuffled);
end
if sum(u==indReport)
    fprintf([num2str(u/10,2),'%%,'])
end
end

if zoom==1

subplot(3,2,4)
    intervall = [-size(shuffled_data,2):10:size(shuffled_data,2)];
histogram = histc(latencies_shuffled,intervall);
bar(intervall,histogram,'k');
xlabel('Latency (ms)')
ylabel('Number of spikes')
xlim([-size(shuffled_data,2) size(shuffled_data,2)])
axis square
box off
    
subplot(3,2,2)
[a,b]=find(shuffled_data==1);
sortiert = sortrows([a,b]);
for j= 1:number
    for i = 1:length(sortiert)
        if sortiert(i,1)==j
            line([sortiert(i,2) sortiert(i,2)],[j-1 j],'Color','k','LineWidth',1)
            xlim([0 length(demean_data_call)/40])
            xlabel('Time (ms)')
            ylabel('Trials')
        end
    end
end
axis square

 subplot(3,2,6)
% % 
 totalspikes= [0:0.1: 100*max(STA_number_shuffled)];

 histogram_totalspikes = histc(100*STA_number_shuffled, totalspikes);
bar(totalspikes,histogram_totalspikes,'k');
 ylabel('number')
 xlabel('% Latencies within +-10ms')
 
 hold on 
 plot(100*STA_number,1,'or','MarkerSize',10,'MarkerFaceColor','r')
box off
 axis square
 
end

totalspikes= [0:0.1: 100*max(STA_number_shuffled)];
sorted_STA_shuffled = sort(STA_number_shuffled);


    viBelow = (sorted_STA_shuffled<=STA_number);
    newSTAShuffled = [sorted_STA_shuffled(viBelow) STA_number sorted_STA_shuffled(~viBelow)];
p_value = 1-max(find(newSTAShuffled==STA_number))/shuffle_time     

zscore_STA=(STA_number-nanmean(STA_number_shuffled))/nanstd(STA_number_shuffled); % calculating precision score
precision_final = sign(zscore_STA)*sqrt(abs(zscore_STA));

zscore_STA_cells(z)=zscore_STA;
precision_final_cells(z)=precision_final;

p_value_all(z)=p_value;
avg_ISI=mean(ISI);
std_ISI=std(ISI);
avg_ISI_silence=mean(ISI_silence);
std_ISI_silence=std(ISI_silence);
STA_number_cells(z)=STA_number; 

elseif length(idx_spikes_col) <= 1
    STA_number= NaN(1);
    STA_number_shuffled=NaN(1);
end
STA_number_cells2(z)=STA_number;

    else
       p_value_all(z)=NaN;
    precision_final_cells(z)=0; 
    end
    else 
    p_value_all(z)=NaN;
    precision_final_cells(z)=0;
    

end



clearvars -except nozerotrials fr_cell_all p_value_all library fields treat zscore_STA_cells precision_final_cells STA_number_cells2 analysis_type timing_type filelist2 number1 fs pre_silence post_silence pathname11 pretime posttime bird bird_ID treatment stdv_fr_call stdv_fr_playback dph cell_ID
end
