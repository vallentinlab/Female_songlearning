
%% Firing rate during female call is measured by duration of female call + 0.05 s:


clear all
close all
clc
%% Load data/ adress folder %%
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

%% choose .mat file that includes aligned sound and recording traces

% [pathname] = uigetdir('DIRECTORY FOR FILES');
% eval(['cd ' pathname]);

% for silence right before call treat is 'sil_on' (if call is before song)

filelist2 = dir('*_E.mat');
treat=('sil_on'); % set treat to 'no' to look at calls only and 'sil_no' to look at silence before

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
avg_fr_call=[];
avg_fr_playback=[];
stdv_fr_call=[];
stdv_fr_playback=[];
dph=[];
cell_ID=[];

%% determine call type and timing:

sil=input('Which silcence: before(1) or after(2) do you want to compare to the playback?');


for z = 1:number1
    
filename = char(strcat(pathname11,'\',filelist2(z,1).name));
load(filename)

%% for juvenile headfixed data:

% bird = (filename(118:121));
% bird_ID = [bird_ID; bird];
% dp = (filename(129:130));
% dph = [dph; dp]
% cell = (filename(123:127));
% cell_ID = [cell_ID; cell];
% old_traces=traces;

%% for juv data:

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


thresh=15;%thresh1(z);

%% detect where song and call is:

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

onset_call_digit = fix(call_onset_ms*40);
offset_call_digit = fix(offset_call_ms*40+2000);%added 50 ms

call_length = fix(offset_call_digit-onset_call_digit); % digits
call_length_s = call_length/40000;
call_length_ms = call_length_s*1000;



playback_matrix = zeros(number, stim_length*40);
call_matrix=zeros(number, call_length*40);
silence_matrix2 = zeros(number,pretime);

positions_silence_before= posttime;
presilence_and_song=length(traces{1,1})-posttime;


for p=1:length(traces)

mVS2=traces{1,p}(:,1)-smooth(traces{1,p}(:,1),fs/80);%high pass filter it. here it is at 2ms 
% thresh = 20; % for juvenile cells
% 
% if z==13
%     thresh = 15; 
% end

g=figure

plot(mVS2)
hold on
plot([1 length(mVS2)], [thresh thresh], '-r')
ylim([-10 70])

pause

close(g)

[pksp,time_spike_playback] = findpeaks(mVS2,'MinPeakHeight',thresh,'MinPeakDistance', 80);
number_of_spikes_playback(p) = length(findpeaks(mVS2,'MinPeakHeight',thresh,'MinPeakDistance', 80));

number_of_spikes_playback2(p) = length(time_spike_playback);

fr_playback(p) = number_of_spikes_playback(p)/stim_length_sec;

end

fr_cell=mean(fr_playback);
fr_cell_all(z)=fr_cell;

clear mVS2 pksp time_spike_playback number_of_spikes_playback fr_playback fr_cell

%% check if celll activity is above 1 Hz

%if fr_cell_all(z)>=1 % check if this cell has at least 1 Hz firing rate

%% compare firing rate control against file with call:

for k = 1:number
    %figure(k+1)
    Fs=40000;

t=(0:length(traces{1,k}(:,1))-1)/Fs;
mVS2=traces{1,k}(:,1)-smooth(traces{1,k}(:,1),Fs/80);%high pass filter it. here it is at 2ms 



if sil == 1
    
    silence=(1:pretime);
    
elseif sil == 2
    
    silence=(pretime+stim2:pretime+stim2+posttime);
    
end

[pksp,time_spike_playback] = findpeaks(mVS2(pretime+1:pretime+1+stim_length*40),'MinPeakHeight',thresh,'MinPeakDistance', 80);
[pksp1,time_spike_call] = findpeaks(mVS2(onset_call_digit:offset_call_digit),'MinPeakHeight',thresh,'MinPeakDistance', 80);
[pksp2,time_spike_silence] = findpeaks(mVS2(silence),'MinPeakHeight',thresh,'MinPeakDistance', 80);


number_of_spikes_playback(k) = length(time_spike_playback);
number_of_spikes_call (k) = length(time_spike_call);
number_of_spikes_silence(k) = length(time_spike_silence);

fr_playback(k) = number_of_spikes_playback(k)/stim_length_sec;
fr_call(k) = number_of_spikes_call(k)/call_length_s;
fr_silence(k) = number_of_spikes_silence(k)/(length(silence)/fs);




end


%% calculate and save values:

all_fr_call{z}=fr_call;
avg_fr_call(z) = mean(fr_call);
stdv_fr_call(z) = std(fr_call);
nr_trials(z) = length(fr_call);
nr_spikes_call{z}=number_of_spikes_call;

sum_spikes_call(z) = sum(number_of_spikes_call);

if strcmp(treat,'no')==1 || strcmp(treat,'sil_no')==1
       
clear thresh bird dp  cell t traces sound stim stim2 stim_length stim_length_sec filename number_of_spikes_playback number_of_spikes_call
clear number call_matrix playback_matrix silence_matrix silence_matrix2 dotraster mVS2 pksp number_of_spikes_silence
clear time_spike_playback time_spike_call number_of_spikes_playback2 number_of_spikes_call2 fr_playback fr_call
clear pks pks_post time_spike_pre_playback time_spike_post_playback number_of_spikes_post_playback number_of_spikes_silence
clear all_spikes all_spikes_psth one_third_stim onset_call_digit call_onset_ms offset_call_digit offset_call_ms
clear positions_silence_before presilence_and_song thresh handle number_of_spikes_call    
    
continue

all_fr_playback{z}=fr_playback;
avg_fr_playback(z) = mean(fr_playback);
stdv_fr_playback(z) = std(fr_playback);
nr_spikes_playback{z}=number_of_spikes_playback;


all_fr_silence{z}=fr_silence;
avg_fr_silence(z) = mean(fr_silence);
stdv_fr_silence(z) = std(fr_silence);
nr_spikes_silence{z}=number_of_spikes_silence;

end

clear thresh bird dp  cell t traces sound stim stim2 stim_length stim_length_sec filename number_of_spikes_playback number_of_spikes_call
clear number call_matrix playback_matrix silence_matrix silence_matrix2 dotraster mVS2 pksp number_of_spikes_silence
clear time_spike_playback time_spike_call number_of_spikes_playback2 number_of_spikes_call2 fr_playback fr_call
clear pks pks_post time_spike_pre_playback time_spike_post_playback number_of_spikes_post_playback number_of_spikes_silence
clear all_spikes all_spikes_psth one_third_stim onset_call_digit call_onset_ms offset_call_digit offset_call_ms
clear positions_silence_before presilence_and_song thresh handle number_of_spikes_call
end
