clear all
close all
clc

%% Load data/ adress folder %%
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

filelist = dir('*.xlsx');
number = length(filelist);
bird_ID = [];
day_nr = [];
treatment = [];
withfem=[];
isolate=[];
c1=0;
c2=0;



for i = 1:number
filename = char(strcat(pathname,'\',filelist(i,1).name));
[data,txt] = xlsread(filename);

bird = (filename(end-23:end-20));
bird_ID = [bird_ID; bird];

treat = (filename(end-11:end-5));
treatment = [treatment; treat];
similarity = data(:,1);
accuracy = data(:,2);

mean_similarity(i)=mean(similarity);

global_similarity(i) = (mean(similarity.*accuracy)/100)./0.7993;
standard_deviation(i) = std(((similarity.*accuracy)/100)./0.7993); 


if strcmpi(treat, 'withfem')
    c1=c1+1;
    plot(1,global_similarity(i),'ro', 'Linewidth',4, 'MarkerSize', 10) 

    fem_idx(c1)=(i); 
    hold on
   
elseif strcmpi (treat,'isolate')
    c2=c2+1;
    plot(2,global_similarity(i), 'bo', 'Linewidth',4, 'MarkerSize', 10)

     iso_idx(c2)=(i);
    hold on
end


end

plot([0.8 1.2],[mean(global_similarity(fem_idx)) mean(global_similarity(fem_idx))], 'k-', 'LineWidth', 4)
hold on
plot([1.8 2.2],[mean(global_similarity(iso_idx)) mean(global_similarity(iso_idx))], 'k-', 'LineWidth', 4)
hold on
xlim([0 4])
axis square
set(gca, 'XTick', [])



box off

ylabel('% Global Similarity')

title('Vocal Learning')
set(gca,'FontSize',30)
set(gcf,'color',[1 1 1])

stats_femiso=ranksum(global_similarity(fem_idx),global_similarity(iso_idx))