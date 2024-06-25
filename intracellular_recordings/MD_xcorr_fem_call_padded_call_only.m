clc
close all
clear all

[pathname2] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname2]);


 load('Combined_alignment_MD_femcall_paddedtraces_09072023.mat')

sma=0;
Call_nr=0;

for p=1:length(padded_traces)



for t=1:length(bird.call_id{p}) % only used to step through all traces
    
    Call_nr=Call_nr+1;
    
    if Call_nr>=16 && Call_nr<21
        
            pretime=100*40; %ms to units
            posttime=50*40;
            
    elseif Call_nr==12
        
            pretime=100*40; %ms to units
            posttime=50*40;
            
    else
        pretime=50*40;
        posttime=50*40;
            
    end
   
    current_traces=padded_traces{t,p};
    current_sound=padded_sound{t,p};
    
    main_trace=current_traces(1,pretime:end-posttime); % call onset to offset
    number=size(current_traces,1);

    
    fem_count=0;
   
    
    %% how many spikes per trace:

    thresh=15;
    
    
    for a=1:number 
        
        fem_count=fem_count+1;
        
        if a==1
            if max(main_trace>thresh)
                
            [pksp_fem,time_spike_fem] = findpeaks(main_trace,'MinPeakHeight',thresh,'MinPeakDistance', 80);
                FEM_spikecount(Call_nr)=length(time_spike_fem);
                
                clear time_spike_fem pksp_fem
            else
                    FEM_spikecount(Call_nr)=0;
            
            end
            
        else
            if max(current_traces(a,pretime:end))>thresh    
            
              [pksp_bos,time_spike_bos] = findpeaks(current_traces(a,pretime:end-posttime),'MinPeakHeight',thresh, 'MinPeakDistance',80);
            BOS_precount{t,p}(a-1)=length(time_spike_bos);
            clear time_spike_bos pksp_bos
            else
                BOS_precount{t,p}(a-1)=0;
            end
        end
        
 
        
    end
    
    BOS_spikecount(Call_nr)=mean(BOS_precount{t,p});
    BOS_spikecount_std(Call_nr)=std(BOS_precount{t,p});
    
    for g=2:number-1
    
        BOS_predeltaspikecount{Call_nr,1}(g-1)=abs(BOS_precount{t,p}(g)-BOS_precount{t,p}(g-1));  
        
    end
    
    mean_predeltaspikecount(Call_nr)=mean(BOS_predeltaspikecount{Call_nr,1});
     std_predeltaspikecount(Call_nr)=std(BOS_predeltaspikecount{Call_nr,1});
 
    Delta_spikecount(Call_nr)=abs(FEM_spikecount(Call_nr)-BOS_spikecount(Call_nr));

    
 figure(3)

 plot(rand, Delta_spikecount(Call_nr),'or')
 hold on
  if Call_nr==22
     
     avg2=mean(mean_predeltaspikecount);
     std2=std(mean_predeltaspikecount);
     
     plot([0 1], [avg2 avg2], '-k')
    hold on 
    plot([0 1], [(1.96*std2)+avg2 (1.96*std2)+avg2], '--k')
    hold on
 end
 xlabel('Call perturbations')
 ylabel('Delta number of spikes')
 axis square
 xlim([-0.5 1.5])
 ylim([-0.2 2.2])
 
 figure(4)
 
 plot(Call_nr, Delta_spikecount(Call_nr),'*b')
 hold on
 xlabel('Call_ID')
 ylabel('Delta abs(FEM_spikes-mean(BOS_spikes))')
 
 if Call_nr==22
     
     avg=mean(Delta_spikecount);
     std=std(Delta_spikecount);
     
     plot([1 22], [avg avg], '-k')
    hold on 
    plot([1 22], [std+avg std+avg], '--k')
    hold on
 end
    


    clear bos_coef b_x_lag b_x_corr bos_vec bos_delta_max b_std b_sem
song_envelope=smooth(abs(hilbert(current_sound(1,:))),50);
    
stch=[1,3,5,7];
atch=[2,4,6,8];

f1=figure

dist_sharp = 30; 

subplot(4,1,1)

    
spectrogram(current_sound(2,:),1024,1000,1024,40000,'yaxis');

colorbar off

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-35 5])
ylim([1 7])
xlim([13 198])
axis square
hold on

if p==10
    caxis([-60 5])
end


for count_r=2:number
    
     subplot(4,1,2)

plot(1/40:1/40:length(current_traces(count_r,:))/40, dist_sharp*(count_r-1)+current_traces(count_r,:), '-', 'Color', [0.7 0.7 0.7])
hold on
plot([pretime/40 pretime/40],[-5 (number-1)*40], '-m', 'LineWidth', 1)
hold on
plot([(length(main_trace)/40)+pretime/40 (length(main_trace)/40)+pretime/40],[-5 (number-1)*40], '-m', 'LineWidth', 1)
hold on
axis square
ylim([20 103])
xlim([13 198])

end

subplot(4,1,3)


spectrogram(current_sound(1,:),1024,1000,1024,40000,'yaxis');

colorbar off

myColorMap = jet; % Make a copy of jet.

% Assign black(all 0's) to black (the first row in myColorMap).
myColorMap(1, :) = [0 0 0];
colormap(myColorMap)
caxis([-35 5])
ylim([1 7])
xlim([13 198])

axis square

if p==10
    caxis([-60 5])
end


subplot(4,1,4)

plot(1/40:1/40:length(current_traces(1,:))/40,current_traces(1,:), '-k')
hold on
plot([pretime/40 pretime/40],[-10 40], '-m', 'LineWidth', 3)
hold on
plot([(length(main_trace)/40)+pretime/40 (length(main_trace)/40)+pretime/40],[-10 40], '-m', 'LineWidth', 3)
ylim([-10 45])
xlim([13 198])

axis square


 end

 close(f1)   

    
    clear current_traces current_sound main_trace number
end


 
figure
plot(BOS_spikecount, 'ok')
hold on
plot(FEM_spikecount, '*r')

