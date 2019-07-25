% Code written by Christian Pedersen
% Michael Bruchas Lab - UW

%% Code purpose

% This code is designed to take a .mat file (from extraction code) and a
% .txt file of time stamps and generate peri-event measurements of average
% photometry activity.

% This code is designed to align photometry traces with specified times of
% behavioral events. Ultimately, a 2D matrix (LickTrig) is generated where
% each row corresponds to fluorescence intensity over a time period surrounding each
% behavioral event. Rows correspond to events, Columns correspond time

%% Reset MatLab workspace - clears all variables and the command window

clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window


%%

k = 1; % session counter

for sess = 1:2  % this is called a "for loop", google this

    
    if sess == 1   % this is called an "if statement", google this
    photoname = 'C89_M1_051418_drug.mat';   % specify .mat file name
    pokeTimes1 = dlmread('89-1_scoretimes.txt'); % specify .txt file name
    end

    if sess == 2
    photoname = 'C89_M2_051418_drug.mat';
    pokeTimes1 = dlmread('89-2_scoretimes.txt');
    end
    

%%  

experDuration = 1800;  % duration of session (seconds)


%%

Fs = 1017.25;

C1 = pokeTimes1;

startdelay = 0;
load(photoname) % load in .mat photometry file

% trim photometry array to length of session
photom1 = data1(round(Fs*startdelay)+1:round(Fs*(experDuration+startdelay)));
time = linspace(1/Fs,experDuration,experDuration*Fs);

behavior = C1;
timeBehav = linspace(1,experDuration,length(behavior));
 
    % Z score photometry data (will override dF/F)
   photom1 = (photom1-mean(photom1))./std(photom1);

for n =1
    
    figure(sess+300)
    
    plot(timeBehav,behavior+10,'g',...
        decimate(time,1000),decimate(photom1,1000));
    
    ylabel('Z score                       Events per second')
    xlabel('Time (sec)')
    yticks([-10 -5 0 5 10 15 20 25 30 35 40])
    yticklabels({'-10','-5','0','5','0','5','10','15','20','25','30'})
    
    windtop = 20;

     axis([0 experDuration -5 14])

end


%% mathy math

behav1 = behavior;

% if 1 second time bin is True, then flag next 20 time bins
for q = 1:(length(behav1)-20)
    
    if behav1(q) == 1
        
        behavior(q+1:q+20) = 1;
    end
    
end

% flag times when behavioral event starts
n = 1;
for p = 2:length(behavior)

    if behavior(p-1) == 0
        if behavior(p) == 1
        
           licks1(n) = timeBehav(p);
           n = n + 1; 
           
        end
    end
    
end


  
  timewindow = 30; % +/- seconds from event, determines how long to average over

%% Lick triggered averaging (of photom trace)

% optional: downsample photometry trace to make files smoother, smaller
decifactor = 1; % how much to downsample (10 = 10 times fewer samples)
photom1 = decimate(photom1,decifactor);
time = decimate(time,decifactor);
 
% update photometry sample rate with new resampled rate
Fs = 1017.25/decifactor;

samplewindow = round(timewindow*Fs); % samples per time period

% remove behavior time stamps that are out of bounds!
licks2 = licks1(licks1<(experDuration-timewindow));
licks = licks2(licks2>timewindow);

% preallocate Matrix for LickTrig
LickTrig = zeros(length(licks),2*samplewindow);

% identify which array indices are adjacent in time to behavior events
for p = 1:length(licks)
    
    trigindex = zeros(1,length(photom1));
    
    [~,startidx] = min(abs(time-(licks(p)-timewindow)));
    trigindex(startidx:(startidx+size(LickTrig,2)-1)) = 1;
    
    LickTrig(p,:) = photom1(trigindex==1);
    
end

% average peri-event fluorescence
photoPerLick = mean(LickTrig,1);

% error bars
sem = std(LickTrig,0,1)./sqrt(size(LickTrig,1)); % sem = std/sqrt(n)

% event triggered average
figure(sess)

hold on
%eb1 = errorbar(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,sem,'Color',[0.2,0.2,0.5]);
x = decimate(linspace(-timewindow,timewindow,length(photoPerLick)),300);
y = decimate(photoPerLick,300);
eb = decimate(sem,300);
lineProps.col{1} = 'blue';
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-2 2]);
%axis([-timewindow timewindow -1.5 1.5])
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Z score (smoothed)')
hold off

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(sess+100)

hold on
imagesc(linspace(-timewindow,timewindow,length(photoPerLick)),1:size(LickTrig,1),LickTrig)
L = line([0 0],[0 length(licks)+1]);
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Bout Number')
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
%colormap(autumn)
xlim([-timewindow timewindow])
ylim([0 length(licks)+1])
hold off

%%


% store 2D matrix from single session into big stack (combine sessions, mice)
eventStack{k} = LickTrig;

k = k + 1;

end


%% combined plot


superStack = cat(1,eventStack{1}); %itchStack{1},itchStack{2},itchStack{3});

for p = 2:length(eventStack)

    superStack = cat(1,superStack,eventStack{p});
    
end


% superStack = peri-event stack for all sessions, combined
LickTrig1 = superStack;


photoPerLick = mean(LickTrig1,1);
%

sem = std(LickTrig1,0,1)./sqrt(size(LickTrig1,1)); % sem = std/sqrt(n)

% event triggered average
figure(61)

hold on
%eb1 = errorbar(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,sem,'Color',[0.2,0.2,0.5]);
x = decimate(linspace(-timewindow,timewindow,length(photoPerLick)),300);
y = decimate(photoPerLick,300);
eb = decimate(sem,300);
lineProps.col{1} = 'blue';
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-2 2]);
%axis([-timewindow timewindow -1.5 1.5])
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Z score (smoothed)')
hold off

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(62)

hold on
imagesc(linspace(-timewindow,timewindow,length(photoPerLick)),1:size(LickTrig1,1),LickTrig1)
L = line([0 0],[0 size(LickTrig1,1)+1]);
set(L,'Color','black')
xlabel('Peri-Event Time (sec)')
ylabel('Bout Number')
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
%colormap(brewermap_view)
colormap(flipud(brewermap([],'YlGnBu')))
xlim([-timewindow timewindow])
ylim([0 size(LickTrig1,1)+1])
hold off



























