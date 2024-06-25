%% If xcorr is NaN, check spike detection threshold for cutting spikes
clear all
close all
clc
%% Load data/ adress folder %% 
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

%% choose .mat file that includes aligned sound and recording traces

filelist2 = dir('*_B.mat');
treat=('G'); % set treat to 'no' to look at calls only and 'sil_no' to look at silence before
                  % and set it to sil_on if looking at juvenile headfixed E
                  % silence before call

number1 = length(filelist2);

%% load in the onset/offset library:
load('Spike_threshold.mat')
load('Solving_percentage_problem_onset_offset_library.mat');


fields=fieldnames(library);

%% variables:

fs = 40000; % Sampling rate


pre_silence=1000;
post_silence=1000;


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

%% determine call type and timing:

    
for z = 1:number1

    

filename = char(strcat(pathname11,'\',filelist2(z,1).name));
load(filename)

%% for juvenile cells:   
bird = (filename(end-27:end-24));
bird_ID = [bird_ID; bird];
cell = (filename(end-10:end-6));
cell_ID = [cell_ID; cell];
old_traces=traces;


stim = (length(traces{1,1})-pretime-post_silence*40)/40;
stim2 = (length(traces{1,1})-pretime-posttime); % stim length in digits
stim_length = fix(stim);
stim_length_sec =  stim_length/1000;

number = size(traces,2);

lil_avg=mean(traces{1,1});
lil_trace=traces{1,1}-lil_avg;

thresh=15;



if z==13
    
    thresh=20; % for awake juveniles
    
end



clear lil_avg lil_trace

for k = 1:number

temp_avg=mean(traces{1,k}(:));
temp_trace=traces{1,k}(:)-temp_avg;


%% Remove spikes (Jagadeesh 1997) %% removes spikes from all traces
spikes=[];
spikes= find(temp_trace>thresh);
for p = 1 : length(spikes)
        if spikes(p)<121 ||spikes(p)> size(traces{1,k}(:,1),1)-500 
        continue
    end
traces{1,k}([spikes(p)-120:spikes(p)+120],1)=NaN;
end

%% interpolating traces

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

clear spikes temp_trace temp_avg
end



%% detect where song and call is:


for h=1:size(fields,1)

    if bird(end-2:end)==fields{h}(end-2:end) % check if the bird matches the library
        
                if strcmp(treat,'sil_on')==1
    
    call_onset_ms=library.(fields{h}).onset.(treat)-100 % take a 50 ms gap between silence before offset and call to avoid fake precision
    offset_call_ms=library.(fields{h}).offset.(treat)-100 % and still have space for +50 ms for auditory response later on
    
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

call_length = offset_call_digit-onset_call_digit; % digits
call_length_s = call_length/40000;
call_length_ms = call_length_s*1000;

%% continue the plot

subplot(4,1,4)

cellmat=cell2mat(traces_new).';
mean_vec=mean(cellmat);
plot(1/40:1/40:length(mean_vec)/40,mean_vec,'b')
hold on
plot([onset_call_digit/40 onset_call_digit/40], [min(mean_vec) max(mean_vec)], '-r')
hold on
plot([offset_call_digit/40 offset_call_digit/40], [min(mean_vec) max(mean_vec)], '-m')
xlim([1/40 length(mean_vec)/40])

subplot(4,1,1)

spectrogram((sound{1,1}),1024,1000,1024,40000,'yaxis');
colorbar off
ylabel({'Frequency', '(Hz)'})% make it in two lines with cell array

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
ylim([1 7])
caxis([-50 0])
sound_envelope=smooth(abs(hilbert(sound{1,1})),50);

clear cellmat 



%% detect where silence is:

silence=1:pretime;



%% De-mean data by substracting the mean value of each trace %% demeans all traces
for r = 1:number
averagemembranepotential_call = mean(traces_new{1,r}(onset_call_digit:offset_call_digit,1),1); %% why minus 1000?
averagemembranepotential_playback=mean(traces_new{1,r}(pretime:pretime+stim2,1),1);
averagemembranepotential_silence=mean(traces_new{1,r}(silence,1),1);


demean_data_call(r,:) = traces_new{1,r}(onset_call_digit:offset_call_digit,1)-averagemembranepotential_call;
demean_data_playback(r,:)=traces_new{1,r}(pretime:pretime+stim2,1)-averagemembranepotential_playback;
demean_data_silence(r,:)=traces_new{1,r}(silence,1)-averagemembranepotential_silence;

traces_new_call{1,r}=demean_data_call(r,:).';

traces_new_playback{1,r}=demean_data_playback(r,:).';

traces_new_silence{1,r}=demean_data_silence(r,:).';
end


clear k



nr_trials(z)=length(traces_new);


%% Calculate xcorr


traces_new_call_mat=cell2mat(traces_new_call).';
traces_new_playback_mat=cell2mat(traces_new_playback).';
traces_new_silence_mat=cell2mat(traces_new_silence).';



coefficient_cal=[];
coefficient_call=[];
coefficient_play=[];
coefficient_silence=[];
for m = 1:number
    
    for n = 1:number
        
        if n == m || n<m
            continue
        end

        crosscorrall_call = xcorr(traces_new_call_mat(m,:),traces_new_call_mat(n,:),0,'coeff'); %40 000 stands for silent period (1s) before and after stimulus

         crosscorrall_play = xcorr(traces_new_playback_mat(m,:),traces_new_playback_mat(n,:),0,'coeff');

          crosscorrall_silence = xcorr(traces_new_silence_mat(m,:),traces_new_silence_mat(n,:),0,'coeff');
         
        correlationcoefficients = crosscorrall_call;
        correlationcoefficients_play = crosscorrall_play;
        correlationcoefficients_silence = crosscorrall_silence;
 
        coefficient_cal = [coefficient_cal correlationcoefficients];
        coefficient_play = [coefficient_play correlationcoefficients_play];
        coefficient_silence = [coefficient_silence correlationcoefficients_silence];
        
    end
end



mean_call_trace{z,1}=mean_vec(onset_call_digit:offset_call_digit);

coefficient_call=  coefficient_cal;
avg_corr_call = mean(coefficient_call);
std_corr_call = std(coefficient_call);

coefficient_call_all{z,1}=coefficient_call;
avg_corr_call_all(z)=avg_corr_call;
std_corr_call_all(z)=std_corr_call;



coefficient_play=  coefficient_play;
avg_corr_play = mean(coefficient_play);
std_corr_play = std(coefficient_play);

coefficient_play_all{z,1}=coefficient_play;
avg_corr_play_all(z)=avg_corr_play;
std_corr_play_all(z)=std_corr_play;



coefficient_silence=  coefficient_silence;
avg_corr_silence = mean(coefficient_silence);
std_corr_silence = std(coefficient_silence);

coefficient_silence_all{z,1}=coefficient_silence;
avg_corr_silence_all(z)=avg_corr_silence;
std_corr_silence_all(z)=std_corr_silence;


figure

plot(mean_call_trace{z,1})

 pause
close all

clearvars -except coefficient_call_all avg_corr_call_all std_corr_call all coefficient_play_all avg_corr_play_all...
    std_corr_play_all coefficient_silence_all avg_corr_silence_all std_corr_silence_all z pathname11 filelist2...
    fields number1 fs pre_silence post_silence sil bird bird_ID treatment treat dp dph cell cell_ID...
    pretime posttime i library thresh1 mean_call_trace

end
