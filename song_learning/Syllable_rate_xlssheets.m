clear all
close all
clc

load('Syllrate_workspace_OfersIsolates.mat')

Ofer_syllrate=mean(syll_rate);
Ofer_sem=1.96*((std(syll_rate))/sqrt(length(syll_rate)));

clearvars -except Ofer_syllrate Ofer_sem

%% Load data/ adress folder %%


[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

filelist = dir('*.xlsx');

tut_syll_rate=5.8724; % normal ABs

tut_sylldur_avg=0.1703;
number = length(filelist);
bird_ID = [];
day_nr = [];
treatment = [];
withfem=[];
isolate=[];
c1=0;
c2=0;

for i = 1:length(filelist)
filename = char(strcat(pathname,'\',filelist(i,1).name));
[data,txt] = xlsread(filename);


treat = (filename(end-11:end-5)); %for normal ABs

bird=(filename(end-8:end-5)); % for bird fathers
bird_ID = [bird_ID; bird];

syll_duration{1,i}=data(:,4)-data(:,3);

for z=2:size(data,1)
   
    gap_duration{1,i}(z-1)=data(z,3)-data(z-1,4)
    
end


% calculate syllable duration per file:

for t=1:length(un)
    
    f=find(file_id{1,i}==un(t));
    
    syll_rate_per_file(t)=(length(syll_duration{1,i}(f))/sum(syll_duration{1,i}(f)));
    
    clear f
    
end

%%
std_bird(i)=std(syll_rate_per_file);

sem_bird(i)=std(syll_rate_per_file)/sqrt(length(syll_rate_per_file));

SyllDur_avg(i)=mean(syll_duration{1,i});
GapDur_avg(i)=mean(gap_duration{1,i});

syll_rate(i)=(length(syll_duration{1,i})/sum(syll_duration{1,i})); %SAP tables give syll_dur in ms so *1000 to make it in s

syll_rate_per_file_bird{1,i}=syll_rate_per_file;


clear data syll_rate_per_file
end

