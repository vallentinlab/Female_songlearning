clear all
close all
clc
%% Load data/ adress folder %%
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

filelist_juv = dir('*chan4.wav');

filelist_csv_juv=dir('*chan4.csv');
filelist_csv_fem=dir('*chan0.csv');

for grab=16%1:length(filelist_csv_juv)

% load juvenile data:

filename_csv_juv = char(strcat(pathname11,'\',filelist_csv_juv(grab,1).name));

[juvenile, juv_labels, juv_raw]=xlsread(filename_csv_juv);

juv_start=juvenile(:,1);
juv_end=juvenile(:,2);
juv_label=juv_labels(2:end,1);

%dph(grab)=str2num(filename_csv_juv(123:124)) %(filename_csv_juv(132:133))
dph(grab)=str2num(filename_csv_juv(111:112))
clearvars -except pathname11 juv_start juv_end juv_label dph grab filelist_csv_juv filelist_juv filelist_csv_fem...
    nrcalls_insong song_percent window_percent nrcalls_inwindow act_nrcalls_insong act_song_dur...
    act_nrcalls_aftersong nr_songs act_nrcalls_inwindow nr_interrupted_songs

% .csv files are organized in sequences

% load female data:

filename_csv_fem = char(strcat(pathname11,'\',filelist_csv_fem(grab,1).name));

[fem, fem_labels, fem_raw]=xlsread(filename_csv_fem);


fem_start=fem(:,1);
fem_end=fem(:,2);
fem_label=fem_labels(2:end,1);

% get the file duration

filename_juv = char(strcat(pathname11,'\',filelist_juv(grab,1).name));

[y,fs]=audioread(filename_juv);

filedur=length(y)/fs;

clearvars -except pathname11 juv_start juv_end juv_label fem_start fem_end fem_label filedur...
    dph grab filelist_csv_juv filelist_juv filelist_csv_fem...
    nrcalls_insong song_percent window_percent nrcalls_inwindow act_nrcalls_insong act_song_dur...
    act_nrcalls_aftersong nr_songs act_nrcalls_inwindow nr_interrupted_songs


%% figure out which files are syllables:

count=0;

for a=1:length(juv_label)
   
    if strcmp(juv_label{a},'syllable')==1 % use for proper labels
        count=count+1;
        syll_idx(count,1)=a;
    end
    
end
clear count

syllable_start=juv_start(syll_idx);
syllable_end=juv_end(syll_idx);

clear juv_start juv_end

%% figure out which are calls:

count=0;

for a=1:length(fem_label)
   
    if strcmp(fem_label{a},'call')==1
        count=count+1;
        call_idx(count,1)=a;
    end
    
end
clear count

call_start=fem_start(call_idx);
call_end=fem_end(call_idx);

clear juv_start juv_end


%% see where song is:

count1=1;
count=0;
song_start(1)=syllable_start(1); % the first syllable is the first start of song

for b=1:length(syll_idx)-1
    
    current_gap=abs(syllable_start(b+1)-syllable_end(b)); % look at absolute value to manage several files in one list
    
    if current_gap>0.350
        count=count+1;
        count1=count1+1;
        song_end(count,1)=syllable_end(b);
        song_start(count1,1)=syllable_start(b+1);
    end
    
end

song_end(count+1,1)=syllable_end(length(syllable_end)); % the last syllable has to be the end of song
song_duration=song_end-song_start;
cut=find(song_duration<0.4); % find where song is shorter than 400 ms
song_duration(cut)=[]; % only save songs that are longer than 400 ms
song_start(cut)=[];
song_end(cut)=[];
sum_songdur=sum(song_duration);
song_percent(grab)=(sum_songdur/filedur)*100;
clear count count1



%% do a song related call count (within the defined window):

calls_in_window=[];

count=0;
count_songs=0;
interrupted_song_starts=0;

for u=1:length(call_idx)
    
    temp_onset=call_start(u);
    
    for z=1:length(song_start) % can exchange for song_window_on
   
        if temp_onset>song_start(z) && temp_onset<song_end(z)+0.35 % add 350 ms
            
        count_songs=count_songs+1;
        
        interrupted_song_starts(count_songs)=song_start(z);
    
            count=count+1;
            calls_in_window(count)=u;
            calls_in_window_songid(count)=z;
            
        end
    end
    clear temp_onset
end

interr=unique(interrupted_song_starts);
nr_interrupted_songs(grab)=length(interr);

%% now introduce if female calls happen between songs:

call_before_song_idx=[];
call_in_song_idx=[];
call_after_song_idx=[];



count=0;
count1=0;
count2=0;
count3=0;



for u=1:length(call_idx)
   
    temp_onset=call_start(u);
    
    for g=1:length(song_start)
        

       
        if temp_onset>song_start(g)&& temp_onset<song_end(g)
            count=count+1;
            call_in_song_idx(count)=u;
            
        elseif temp_onset>song_end(g)&& temp_onset<song_end(g)+0.35 % check if there are calls 350 ms after song end
            count2=count2+1;
            call_after_song_idx(count2)=u;
            
        elseif temp_onset<song_start(g) && temp_onset> (song_start(g)-0.35) % check if there are calls 350 ms before song
            
            count3=count3+1;
            call_before_song_idx(count3)=u;
            
        end
    end
    
    clear temp_onset
end



% substract the calls in song from all calls:

calls_outside_song=1:length(call_start);
calls_outside_song(call_in_song_idx)=[];

call_durations=call_end-call_start;

sum_callsbeforesong_dur=sum(call_durations(call_before_song_idx));
sum_callsinsong_dur=sum(call_durations(call_in_song_idx));
sum_callsaftersong_dur=sum(call_durations(call_after_song_idx));

aftersong_dur=0.35*length(song_duration);
beforesong_dur=0.35*length(song_duration);

clear count count1




act_nrcalls_insong(grab)=length(call_in_song_idx);
act_nrcalls_aftersong(grab)=length(call_after_song_idx);
act_song_dur(grab)=sum_songdur;
nrcalls_insong(grab)=(length(call_in_song_idx)/length(call_idx))*100;
nrcalls_inwindow(grab)=(length(calls_in_window)/length(call_idx))*100;
nr_songs(grab)=length(song_start);
act_nrcalls_inwindow(grab)=length(calls_in_window);


clearvars -except window_percent song_percent nrcalls_insong pathname11 dph filelist_juv...
    filelist_csv_juv filelist_csv_fem grab nrcalls_inwindow act_nrcalls_insong act_song_dur...
    act_nrcalls_aftersong nr_songs act_nrcalls_inwindow nr_interrupted_songs

end
