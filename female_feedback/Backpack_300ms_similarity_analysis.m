clear all
close all
clc

%% Load data/ adress folder %%
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

filelist = dir('*.xlsx');
number = length(filelist);
day_nr = [];
treatment = [];

count1=0;
count2=0;
int_sim=[];
unt_sim=[];
nr_snips=0;

for i = 1:number
filename = char(strcat(pathname,'\',filelist(i,1).name));
[data,txt] = xlsread(filename);
%bird = (filename(87:90)); % for backpacks
day_nr = (filename(96:97)); % for non-contingent femstack server
dph(i) = str2num(day_nr);
treatment=(filename(102:112));

similarity = data(:,1);
accuracy = data(:,2);



local_similarity{1,i}=((similarity.*accuracy)/100)./0.7993;

global_similarity(i) = (mean(similarity.*accuracy)/100)./0.7993;

median_birds(i) = (median(similarity.*accuracy)/100)./0.7993;

nr_snips=nr_snips+length(local_similarity{1,i});

std_sim(i)=std(((similarity.*accuracy)/100)./0.7993);

sem(i)=std(((similarity.*accuracy)/100)./0.7993)/sqrt(length(similarity));
all_day = str2num(day_nr);

cat=repmat(dph(i),1,length(local_similarity{1,i}));
ra=rand(1,length(local_similarity{1,i}));

if strcmpi(treatment, 'interrupted')
    count1=count1+1;
    int_idx(count1)=i;
    int_sim=vertcat(int_sim, local_similarity{1,i});
    plot(cat+rand+1,local_similarity{1,i},'ro', 'Linewidth',2)
    hold on
    plot([dph(i)+1 dph(i)+2], [global_similarity(i) global_similarity(i)], '-k', 'Linewidth',2)
    hold on

elseif strcmpi (treatment,'uninterrupt')
    count2=count2+1;
    unt_idx(count2)=i;
    unt_sim=vertcat(unt_sim, local_similarity{1,i});
    plot(cat-rand-1,local_similarity{1,i}, 'bo', 'Linewidth',2)
    hold on
    plot([dph(i)-1 dph(i)-2], [global_similarity(i) global_similarity(i)], '-k', 'Linewidth',2)
    hold on

end



clearvars -except local_similarity global_similarity filelist number pathname unt_sim int_sim...
    int_idx unt_idx count1 count2 dph nr_snips
end
xlabel('dph')
xticklabels(dph)
ylabel('Similarity to tutor song (%)')
axis square
title('Bird 23184')
ylim([0 40])


figure

h1=histogram(int_sim, 10)
hold on
h2=histogram(unt_sim, 10)
ylim([0 10])
xlim([0 40])
xlabel('Similarity to tutor song')
axis square
title('Bird 23184')

alldays_stats=ranksum(int_sim,unt_sim)


for t=1:length(filelist)/2
   
    stats(t)=ranksum(local_similarity{1,int_idx(t)}, local_similarity{1,unt_idx(t)});
    
end