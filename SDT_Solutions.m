% "Solutions" for Day 3 "Analysis of extracellular activity across layers 
% during a perceptual task" of the G-Node Winter course in Neural Data"
% for the GNode Workshop Summer 2017
% This script uses cells that can be run separately. 
% by Hendrikje Nienborg, CIN, Uni Tuebingen
% July, 2017
%
% part 1 (introduction to signal detection theory) is based on Geoffrey
% Boynton's SDT tutorial for the cold spring harbor course
%
%
%% Part 1) -------  Introduction to SDT ----------------------------------
%
%
%% Tasks:
%%  1) variability in the internal response could reflect noise at the 
% level of the sensory representation (e.g. synaptic transmission), 
% or at a higher level, e.g. variability in the subject's motivational, 
% attentional state, behavioral bias, etc.
%
%% 2) Plotting the distributions of the internal response for the noise and 
%% the signal

% the distributions have an SD of 1, and are centered on 0 (for the noise)
% and 1 (for the signal)
noiseMean = 0;   
signalMean = 1;
sd = 1;

% We will calculate the height of each probability distribution at these
% points
z = -4:.2:6;  

noise_y  = normpdf(z,noiseMean,sd);
signal_y = normpdf(z,signalMean,sd);

% plotting the distributions
fh = figure;
subplot(2,1,1)
plot(z,noise_y);
hold on
plot(z,signal_y,'r-');

ylim = get(gca,'YLim');

text(noiseMean,ylim(2)*.9,'Noise','VerticalAlignment','top','HorizontalAlignment','center','Color','b');
text(signalMean,ylim(2)*.9,'Signal','VerticalAlignment','top','HorizontalAlignment','center','Color','r');
xlabel('Internal Response');
ylabel('Probability')

%% 3) exploring the role of the criterion
%
%% a) superimposing the criterion on the plot
criterion = 1; 

plot(criterion*[1,1],ylim,'k:');
hcrit = text(criterion,0,'criterion','VerticalAlignment','bottom','HorizontalAlignment','left');

%% b) in SDT the probabilities for Hits, FAs, CRs, Misses correspond to
% areas under the probability density functions of the internal responses
% for the signal and noise:
%
% pHit: the probability for a Hit corresponds to the area under the
%    signal distribution that exceeds the criterion.
% pMiss: the probability for a Miss corresponds to the area under the
%    signal distribution that is < than the criterion.
% pFA: the probabilities for a Hit corresponds to the area under the
%    noisel distribution that exceeds the criterion.
% pRC: the probabilities for a Hit corresponds to the area under the
%    signal distribution that is < than the criterion.
%%
%% c) computing the probabilities for Hits, FAs, Misses, CRs
pHit  =  1-normcdf(criterion,signalMean,sd)
pFA   =  1-normcdf(criterion,noiseMean,sd)
pMiss =  normcdf(criterion,signalMean,sd)
pCR   =  normcdf(criterion,noiseMean,sd)


%% d1) the criterion should be higher such that pHit>2*pFA
criteria = [0:0.1:3];
pHits = 1-normcdf(criteria,signalMean,sd);
pFAs   =  1-normcdf(criteria,noiseMean,sd);
fh2 = figure;
xlim = [min(criteria) max(criteria)];
plot(criteria,(pHits-(2*pFAs))*10,'bo'); hold on
plot(xlim,zeros(1,2),'--k')
xlabel criterion
ylabel('average win / trial (cents)')

% d2) in the situation for cancer detection the cost of a miss is very
% high. So the criterion should be very low, i.e. there will many FAs.

%% e) instead of the wins we will now maximize the proportion correct. Since
% on half of the trials there is a signal, and on half there is none, we
% will average between pHit and pCR, and check at which criterion this
% value is maximal.
figure(fh)
subplot(2,1,2)
criteria = -4:.1:6; 
pHits = 1-normcdf(criteria,signalMean,sd);
pCRs   =  normcdf(criteria,noiseMean,sd);
proportion_correct = (pHits+pCRs)/2;
plot(criteria,proportion_correct,'bo'); hold on

xlabel criterion
ylabel('proportion correct')

plot([.5 .5] ,get(gca,'ylim'),'--k')

%% 4) The receiver operating characteristic (ROC) curve
%% a) the higher the Hit rate, the higher the FA rate, too.

%% b) plotting the ROC curve for the initial signal and noise distribution
criteria = -4:.1:6; 
pHits = 1-normcdf(criteria,signalMean,sd);
pFAs   =  1-normcdf(criteria,noiseMean,sd);
fh3 = figure;
plot(pFAs,pHits,'-k')
xlabel('FA rate')
hold on
plot([0,1],[0,1],'k:');ylabel('Hit rate')

title (' ROC curve')

%% c) plotting the ROC curve for identical distibutions
signalMean = 0;
noiseMean = 0;
pHits = 1-normcdf(criteria,signalMean,sd);
pFAs   =  1-normcdf(criteria,noiseMean,sd);
fh4 = figure;
plot(pFAs,pHits,'-k')
xlabel('FA rate')
hold on
plot([0,1],[0,1],'k:');ylabel('Hit rate')

title ('ROC curve')

%% d) the aROC
signalMean = 1;
noiseMean = 0;
pHits = 1-normcdf(criteria,signalMean,sd);
pFAs   =  1-normcdf(criteria,noiseMean,sd);
aROC = -trapz(pFAs,pHits);  % -trapz  to account for decreasing FA rates for increasing criteria)

%% 5) simulating a 2AFC experiment. We will show that aROC corresponds to 
%% percent correct in a 2AFC paradigm.
signalMean = 1;
noiseMean = 0;
sd = 1;

nTrials = 10000;

x= randn(2,nTrials); % internal responses for signal and noise intervals
x(1,:) = x(1,:)*sd + noiseMean;
x(2,:) = x(2,:)*sd + signalMean;

response = x(2,:)>x(1,:);

% response is 1 for correct trials, when the draw from the signal exceeds
% the noise draw.
%
% Overall performance is the mean of the response vector:
percentCorrect = mean(response);
fprintf('Percent correct: %5.3f\n',percentCorrect);
fprintf('aROC: %5.3f\n',aROC);

%% 6) optional 1: the distribution for the signal would be to the left of 
%% that for the noise, i.e. signalMean < noiseMean.
%%  optional 2: perform examples in perfcurve


%% Part 2) -- using SDT to analyze neuronal data--------------------------
%% Tasks
%% 1) looking at the data
%% a/b) plotting the histogram and mean +/-SE responses for each stimulus
curr_dir = cd('C:\Users\franz\Documents\gnode\GNode2017\DataGNode');
load ch19.mat
cd(curr_dir);

s = m(:,2); %stimulus
S = unique(s); % stimulus types
xlim = [min(m(:,3)), max(m(:,3))];
bins = [xlim(1):diff(xlim)/20:xlim(2)];
fh5 = figure;
sHi = 1/(length(S)+3);

for n = 1:length(S)
    subplot('position',[.3,1-sHi*(n+.5),0.6,sHi]);
    idx = find(s==S(n)); % indices for trials with this stimulus
    hist(m(idx,3),bins);
    if n<length(S)
        set(gca,'xlim',xlim,'box','off','xticklabel','')
    else
        set(gca,'xlim',xlim,'box','off')
        xlabel ('spike rate')
        ylabel ('N')
    end
    text(xlim(1)-diff(xlim)/3,0,sprintf('Stim: %5.2f\n',S(n)));
    meanRate(n) = mean(m(idx,3));
    SE(n) = std(m(idx,3))/sqrt(length(idx));
end

fh6 = figure;
errorbar(S,meanRate,SE,'ks-');
xlabel('stimulus')
ylabel('mean rate (spikes/sec)')
set(gca,'xlim',[-1.1 1.1])

%% 2) the neurometric function
%% a/b)
N = 100; 
noiseD = m(s==0,3);
fa=[];
hit=[];
aROC_neurometric = [];

figure
subplot(2,1,1)
hist(noiseD)
subplot(2,1,2)
hist(signalD)


    signalD = m(s==.5,3);

    criteria = linspace(min([noiseD;signalD])-1,max([noiseD;signalD]),N);
    for j = 1:N
      fa(j) = sum(noiseD > criteria(j))/length(noiseD);
      hit(j) = sum(signalD > criteria(j))/length(signalD);
    end
    
    figure
    plot(hit,fa)
    xlim([-1 2])
    aROC_neurometric(n) = -trapz(fa,hit);


figure;
subplot(1,2,1)
plot(S,aROC_neurometric,'-bo');
xlabel('stimulus')
ylabel('aROC')
title ('neurometric function')


N = 100; 
noiseD = m(s==0,3);
fa=[];
hit=[];
aROC_neurometric = [];

for n = 1:length(S)
    signalD = m(s==(S(n)),3);

    criteria = linspace(min([noiseD;signalD])-1,max([noiseD;signalD]),N);
    for j = 1:N
      fa(j) = sum(noiseD > criteria(j))/length(noiseD);
      hit(j) = sum(signalD > criteria(j))/length(signalD);
    end
    
    aROC_neurometric(n) = -trapz(fa,hit);
end

fh7 = figure;
subplot(1,2,1)
plot(S,aROC_neurometric,'-bo');
xlabel('stimulus')
ylabel('aROC')
title ('neurometric function')

%% 3) the psychophysical function
propNearChoice = []; % proportion near choice

for n = 1:length(S)
    propNearChoice(n) = sum(m(s==S(n),4)==0)/sum(s==S(n));
    Nsamples(n) = sum(s==S(n)); 
end
subplot(1,2,2)
plot(S,propNearChoice,'-rs')
xlabel ('stimulus')
ylabel('proportion near choice')
title ('psychometric function')

%% 4) choice probabilities


for n = 1:length(S)
    % based on figure 5 we know that we have a near-preferring neuron
    idx = find(s==S(n));
    near_choices = find(m(idx,4)==0);
    far_choices = find(m(idx,4)==1);
    
    signalD = m(idx(near_choices),3); % we define responses on near choices as signal distribution
    noiseD = m(idx(far_choices),3); % we define responses on far choices as noise distribution
    criteria = linspace(min([noiseD;signalD])-1,max([noiseD;signalD]),N);
    for j = 1:N
      fa(j) = sum(noiseD > criteria(j))/length(noiseD);
      hit(j) = sum(signalD > criteria(j))/length(signalD);
    end
    faX(n,:) = fa;
    hitX(n,:) =hit;
    
    CP(n) = -trapz(fa,hit);
end

figure
for i = 1:size(faX,1)
    plot(faX(i,:),hitX(i,:))
    hold on
    
end

fh8 = figure;
plot(S,CP,'ks')
xlabel('stimulus')
ylabel('choice probability')

%% Part 3) -- optional: testing Pitkow read-out and datachallenge
%% 
%% Tasks 
%% 1) 
cc_measured = pi*(CP-0.5)/sqrt(2);  % choice correlation according to eq. 9 in Pitkow

%% 2) let's fit the psychophysical and neurometric functions with a
%% cumulative Gaussian via max likelihood estimation
addpath 'C:\Users\franz\Documents\gnode\GNode2017\helpers'

% we will need the number of samples for each stimulus
for n = 1:length(S)
    Nsamples(n) = sum(s==S(n)); 
end

% to fit the cumulative Gaussians the curves needs to increase from left to
% right
aROC_neurometric_signed = fliplr(aROC_neurometric); 
psychometric = fliplr(propNearChoice);
Nsamples_signed = fliplr(Nsamples);

x = S'; % make sure, x, y, n have the same size
n = Nsamples_signed;
y = aROC_neurometric_signed;

[fitNeurometric] = FitCumGaussMaxLi(x,y.*n,n,'basefix',0,'ampfix',1); % we force the fit to go from 0 to 1

y = psychometric;
[fitPsychometric] = FitCumGaussMaxLi(x,y.*n,n,'basefix',0,'ampfix',1);% we force the fit to go from 0 to 1

% visualizing the fits
fh9 = figure;
xvals = [-1:.1:1];

subplot(1,2,1)
plot(x,aROC_neurometric_signed,'bo');
hold on; 
plot(xvals,normcdf(xvals,fitNeurometric.mean,fitNeurometric.sd),'-b');
xlabel('stimulus (signed by neural preference)')
ylabel('proportion preferred choice')
title('neurometric')

subplot(1,2,2)
plot(x,psychometric,'ro');
hold on; 
plot(xvals,normcdf(xvals,fitPsychometric.mean,fitNeurometric.sd),'-r');
xlabel('stimulus (signed by neural preference)')
ylabel('proportion preferred choice')
title('psychophysical')

cc_opt = fitPsychometric.sd./fitNeurometric.sd; % predicted optimal choice correlation
fprintf('measured choice correlation: %5.3f\n',cc_measured(7));  % restrict to cc for stimulus = 0 to have balanced choices
fprintf('predicted optimal choice correlation: %5.3f\n',cc_opt);

%% 3) do the above analysis for all the units
%% since all the recordings come from the same session the behavior and the 
%% stimulus remain the same and we only have to re-do the computation for the
%% neuronal data
fh10 = figure;
channels = [1:24];
for ch = channels
    curr_dir = cd('C:\Users\franz\Documents\gnode\GNode2017\DataGNode');
    load(['ch' num2str(ch)]);
    cd(curr_dir);
    % compute the neurometric funtion
    N = 100; 
    noiseD = m(s==0,3);
    fa=[];
    hit=[];
    aROC_neurometric = [];

    for n = 1:length(S)
        signalD = m(s==(S(n)),3);

        criteria = linspace(min([noiseD;signalD])-1,max([noiseD;signalD]),N);
        for j = 1:N
          fa(j) = sum(noiseD > criteria(j))/length(noiseD);
          hit(j) = sum(signalD > criteria(j))/length(signalD);
        end

        aROC_neurometric(n) = -trapz(fa,hit);
    end

    % do we have a "near-preferring" or "far-preferring" unit?
    % we determine the neuronal preference based on the responses to the
    % strongest signal strength (near or far signal)
    if mean(m(s>=0.5,3)) <mean(m(s<=-0.5,3))
        near_flag = 1;
        aROC_neurometric_signed = fliplr(aROC_neurometric);
        n = fliplr(Nsamples);
    else
        near_flag = 0;
        aROC_neurometric_signed = (aROC_neurometric);
        n = Nsamples;
    end
    
    x = S'; % make sure, x, y, n have the same size
    y = aROC_neurometric_signed;

    [fitNeurometric,~,exitflag] = FitCumGaussMaxLi(x,y.*n,n,'basefix',0,'ampfix',1); % we force the fit to go from 0 to 1
    
    % visualize the data and the fits
    subplot(4,6,ch)
    plot(x,y,'bo');
    hold on
    plot(xvals,normcdf(xvals,fitNeurometric.mean,fitNeurometric.sd),'-b');
    title(sprintf('channel: %d', ch))
    
    if ~exitflag
        fitNeurometric.sd = NaN;
    end
    
    % calculate choice probabilities
    % we restrict our calculation of CP to the 0 stimulus
    N=100;
    idx = find(s==0);
    near_choices = find(m(idx,4)==0);
    far_choices = find(m(idx,4)==1);
    if near_flag
        signalD = m(idx(near_choices),3); % we define responses on near choices as signal distribution
        noiseD = m(idx(far_choices),3); % we define responses on far choices as noise distribution
    else
        signalD = m(idx(far_choices),3); % we define responses on far choices as signal distribution
        noiseD = m(idx(near_choices),3); % we define responses on near choices as noise distribution
    end
    
    criteria = linspace(min([noiseD;signalD])-1,max([noiseD;signalD]),N);
    for j = 1:N
      fa(j) = sum(noiseD > criteria(j))/length(noiseD);
      hit(j) = sum(signalD > criteria(j))/length(signalD);
    end
    
    CP = -trapz(fa,hit);
    
    % store the values in a vector
    neuro_sd(ch) = fitNeurometric.sd;
    cc_measured(ch) = pi*(CP-0.5)/sqrt(2); % measured choice correlation
    cc_opt(ch) = fitPsychometric.sd./fitNeurometric.sd;  % predicted choice correlation for optimal linear read-out

end

fh11 = figure;
plot(cc_measured,cc_opt,'ko');
hold on;
unity
xlabel('measured choice correlation')
ylabel('predicted choice correlation')

%% Part 4) -- CSD maps
%% 1) E.g. plotting as colorplot but sorted by channel
% load the lfp datafile
% % it is a matrix lfp = [c by t by n] of the voltage signals 
%where  c: channel #, 
%       t: time samples, [200ms before to 200ms after stimulus onset]
%       n: number of trials, aligned on trial onset, 
% sampling rate is 1kHz

curr_dir = cd('C:\Users\franz\Documents\gnode\GNode2017\DataGNode');
load(['lfp']);
cd(curr_dir)
%%
% let's visualize the first 100 trials for each channel
lfp_100 = lfp(:,:,1:100);
fh12 = figure;

for n = 1:24
    lfp_reshape((n-1)*100+1:n*100,:) = squeeze(lfp_100(n,:,:))';
end

% plotting
subplot(1,5,1)
imagesc(lfp_reshape)
set(gca,'ytick',[50:100:2350],'yticklabel',[1:24],...
    'xtick',[100:100:300],'xticklabel',[-100:100:100])
grid on
xlabel('time rel. to stimulus onset [ms]')
ylabel('channel number')

%% 2) compute the average stimulus-triggered lfp on each channel
lfp_mn = squeeze(mean(lfp,3)); % average stimulus aligned lfp on each channel

subplot(1,5,2)
imagesc(lfp_mn);
set(gca,'ytick',[1:24],'yticklabel',[1:24],...
    'xtick',[100:100:300],'xticklabel',[-100:100:100])
grid on
xlabel('time rel. to stimulus onset [ms]')
title('stimulus triggered lfp per channel')

%% 3) compute 'raw' csd map

csd_raw =  -diff(lfp_mn,2,1); % we compute 2nd spatial derivative

% plotting results
subplot(1,5,3)
imagesc(csd_raw);
set(gca,'ytick',[1:24],'yticklabel',[2:23],...
    'xtick',[100:100:300],'xticklabel',[-100:100:100])
grid on
xlabel('time rel. to stimulus onset [ms]')
title('raw csd map')


%% 4) compute temporally, spatially filtered CSD map, and apply Vaknin method
% compute the default filters
sF=1000; % sampling frequency of our signal

% temporal filters
[low.num low.denom] = butter(4,100/(sF/2),'low'); % 
[high.num high.denom] = butter(2,8/(sF/2),'high');
% spatial filter
spat_filter = fspecial('gaussian',[3 5],1.1);

%filter temporally
lfptmp = filter(high.num, high.denom,lfp_mn,[],2);  
lfptmp = filter(low.num, low.denom,lfptmp,[],2);

% filter spatially
lfptmp = [lfptmp(1,:);lfptmp;lfptmp(end,:)];
lfp = conv2(lfptmp,spat_filter,'same');

% compute the CSD map
lfp_vaknin = [lfp(1,:);lfp;lfp(end,:)]; % vaknin transform
csd = -diff(lfp_vaknin,2,1);

% replace first and last channels of csd map with 0 to ensure same size as lfp
csd = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))]; 
lfp = lfp(2:end-1,:);

% plotting the data
subplot(1,5,4)
imagesc(lfp)
set(gca,'ytick',[1:24],'yticklabel',[1:24],...
    'xtick',[100:100:300],'xticklabel',[-100:100:100])
title('filtered lfp')
grid on
xlabel('time rel. to stimulus onset [ms]')

subplot(1,5,5)
imagesc(csd);
set(gca,'ytick',[1:24],'yticklabel',[1:24],...
    'xtick',[100:100:300],'xticklabel',[-100:100:100])
grid on
xlabel('time rel. to stimulus onset [ms]')
title('filtered, vaknin transformed csd map')
