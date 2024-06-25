%% If xcorr is NaN, check spike detection threshold for cutting spikes
clear all
close all
clc
%% Load data/ adress folder %% 
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

%% choose .mat file that includes aligned sound and recording traces

filelist2 = dir('*b.mat');

number1 = length(filelist2);


%% variables:

fs = 40000; % Sampling rate

%silence before BOS:

pre_silence=500;
post_silence=500;

pretime = pre_silence*40;
posttime=post_silence*40;
i= 1;

bird=[];
bird_ID = [];
treatment = [];
stdv_fr_call=[];
stdv_fr_playback=[];
bird_dph=[];
dph=[];


%% determine call type and timing:

    
for z = 1:number1

    

filename = char(strcat(pathname11,'\',filelist2(z,1).name));
load(filename)

%% for MD juvenile cells:   
bird = (filename(end-32:end-27));
bird_ID = [bird_ID; bird];
dph = (filename(end-25:end-24));
bird_dph = [bird_dph; dph];
old_traces=traces;

number = size(traces,2);

thresh=15; % threshold for spike detection

if z==5
   
    traces_mean=mean(traces{1,1});
    insert(1:pretime,1)=traces_mean;

    
    for u=1:number
       
       traces_new{1,u}=[insert; traces{1,u}(660*40:840*40); insert]; 
        
    end

    clear traces
    traces=traces_new;
    
end


for k = 1:number

temp_avg=mean(traces{1,k}(:));
temp_trace=traces{1,k}(:)-temp_avg; % demean the traces here to cut spikes


%% Remove spikes (Jagadeesh 1997) %% removes spikes from all traces
spikes=[];
spikes= find(temp_trace>thresh);
for p = 1 : length(spikes)
        if spikes(p)<121 ||spikes(p)> size(traces{1,k}(:,1),1)-500  
        continue
    end
traces{1,k}([spikes(p)-120:spikes(p)+120],1)=NaN;
end


% interpolating traces

if isempty(spikes)==0 %check if the spikes were detected and cut
  traces_new{1,k}(:,1)=fixgaps(traces{1,k}(:,1));
  
else
    traces_new{1,k}(:,1)=traces{1,k}(:,1);
end  


figure(z)

subplot(4,1,2)

plot(1/40:1/40:size(temp_trace,1)/40,temp_trace,'m')
xlim([1/40 size(temp_trace,1)/40])

subplot(4,1,3)

plot(1/40:1/40:size(traces_new{1,k},1)/40,traces_new{1,k},'k')
hold on
xlim([1/40 size(temp_trace,1)/40])


end

clear spikes temp_trace temp_avg

%% detect where song and call is:

song_onset_ms=pre_silence+1;
song_offset_ms=(length(traces{1,1})-posttime)/40;


onset_song_digit = song_onset_ms*40;
offset_song_digit = song_offset_ms*40;

song_length = offset_song_digit-onset_song_digit; % digits
song_length_s = song_length/40000;
song_length_ms = song_length_s*1000;

%% continue the plot


subplot(4,1,4)

cellmat=cell2mat(traces_new).';
mean_vec=mean(cellmat);
plot(1/40:1/40:length(mean_vec)/40,mean_vec,'b')
hold on
plot([onset_song_digit/40 onset_song_digit/40], [min(mean_vec) max(mean_vec)], '-r')
hold on
plot([offset_song_digit/40 offset_song_digit/40], [min(mean_vec) max(mean_vec)], '-m')
xlim([1/40 length(mean_vec)/40])

ram=figure(30)

for t=1:number

subplot(number*2,1,t)

spectrogram((sound{1,t}),1024,1000,1024,40000,'yaxis');
colorbar off
ylabel({'Frequency', '(Hz)'})% make it in two lines with cell array

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
ylim([1 7])
caxis([-50 0])
sound_envelope=smooth(abs(hilbert(sound{1,t})),50);

clear cellmat 

subplot(number*2,1,number+t)

plot(1/40:1/40:size(traces_new{1,t},1)/40,traces_new{1,t},'k')
hold on
plot([onset_song_digit/40 onset_song_digit/40], [min(mean_vec) max(mean_vec)], '-r')
hold on
plot([offset_song_digit/40 offset_song_digit/40], [min(mean_vec) max(mean_vec)], '-m')
xlim([1/40 length(mean_vec)/40])
xlim([1/40 size(traces_new{1,t},1)/40])

end

pause
close(ram)
%% detect where silence is:

silence=1:pretime;

%% De-mean data by substracting the mean value of each trace %% demeans all traces
for r = 1:number
averagemembranepotential_song = mean(traces_new{1,r}(onset_song_digit:offset_song_digit,1),1);


demean_data_song(r,:) = traces_new{1,r}(onset_song_digit:offset_song_digit,1)-averagemembranepotential_song;

traces_new_song{1,r}=demean_data_song(r,:).';

end


clear k


%% Calculate xcorr

traces_new_song_mat=cell2mat(traces_new_song).';



coefficient_son=[];
coefficient_song=[];

for m = 1:number
    
    for n = 1:number
        
        if n == m || n<m
            continue
        end

        crosscorrall_song = xcorr(traces_new_song_mat(m,:),traces_new_song_mat(n,:),0,'coeff'); 

        correlationcoefficients = crosscorrall_song;
        coefficient_son = [coefficient_son correlationcoefficients];

    end
end



mean_song_trace{z,1}=mean_vec(onset_song_digit:offset_song_digit);

coefficient_song=  coefficient_son;
avg_corr_song = mean(coefficient_song);
std_corr_song = std(coefficient_song);

coefficient_song_all{z,1}=coefficient_song;
avg_corr_song_all(z)=avg_corr_song;
std_corr_song_all(z)=std_corr_song;




clearvars -except coefficient_call_all avg_corr_call_all std_corr_call all coefficient_play_all avg_corr_play_all...
    std_corr_play_all coefficient_silence_all avg_corr_silence_all std_corr_silence_all z pathname11 filelist2...
    fields number1 fs pre_silence post_silence sil bird bird_ID treatment treat dp dph cell cell_ID...
    pretime posttime i library thresh1 mean_call_trace dph bird_dph coefficient_song_all avg_corr_song_all...
    std_corr_song_all

end

figure
a=repmat(1,9,1);

scatter(a,avg_corr_song_all,'mo','Jitter','on')
hold on
plot([0.8 1.2], [mean(avg_corr_song_all) mean(avg_corr_song_all)], '-k')
ylim([0 1])
xlim([0 2])
axis square
box off
