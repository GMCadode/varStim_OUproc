% script for generating jittered timing for stim pulses
% code for generating stim file of variable IPI stim pulses for Lankarany
% collab


%% Generate entire stim file with four components:
% 15 sec 130 Hz jittered
% 15 sec 130 Hz isochronous
% 15 sec 185 Hz jittered
% 15 sec 185 Hz isochronous
%  
% "jittered" pulse events according to OU process


totaltime = 15; % seconds; total time when stims will be on
    srate = 24414.0625; % samples/sec; system sampling rate of TDT setup (e.g. 24414.0625)
percStdms = 0.1; % scalar; the percentage around avISIms for 1 standard deviation (0.1),
% 10% of a milliseconds

% 
stim = [];

% create first stim file based on jittered 130Hz stim
dbsfreq = 130; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stim_130j = genJitterStimtrain(totaltime, srate, avISIms, percStdms);

% create first stim file based on isochronous 130Hz stim
dbsfreq = 130; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stim_130r = genRegularStimtrain(totaltime, srate, avISIms, percStdms);

% create first stim file based on isochronous 130Hz stim
dbsfreq = 185; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stim_185j = genJitterStimtrain(totaltime, srate, avISIms, percStdms);

% create first stim file based on isochronous 130Hz stim
dbsfreq = 185; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stim_185r = genRegularStimtrain(totaltime, srate, avISIms, percStdms);

stim = [stim_130j; stim_130r; stim_185j; stim_185r];
% t = (0:(length(stim)-1)) * (1/srate);
% figure; plot(t,stim)



% convolve the events with the stimulus waveform



% save the results as a .mat file to be loaded into TDT Synapse software


%%
% fs = 24414.0625;
% 
% % initialize a vector of values that comprise the stimulation file; default
% % zero values
% totRuntime = 60; % seconds
% nSamples = round(fs * totRuntime);
% stimEvents = zeros(nSamples,1);
% 
% % generate a series of pulse time occurrences based on a Poisson process
% avStimRate = 100; % Hz, i.e. pulses per second
% dt = 1/avStimRate
% dS = round(fs * dt);
% lambda = dS; % lambda here pertains to the average stim rate, say 100 Hz
% r = poissrnd(lambda, 12000, 1);
% 
% % space out all the 1's according to the poisson process sample nums
% rCum = cumsum(r);
% rCum(rCum > nSamples) = [];
% stimEvents(rCum) = 1;


% as a safeguard measure, control for all ISI's that would cause the pulse
% waveform to continue on to the next one too


% Convolve stim pulse waveform with event occurrences
stimFile = conv(stimEvents, waveform); 
lWav = length(waveform);
stimFile(end-(lWav-2):end) = [];
tFile = ((1:nSamples)-1) * (1/fs);
figure; plot(tFile, stimFile)



%% Old code
% 
% % script for generating jittered timing for stim pulses
% 
% fs = 24414.0625;
% 
% % initialize a vector of values that comprise the stimulation file; default
% % zero values
% totRuntime = 60; % seconds
% nSamples = round(fs * totRuntime);
% stimEvents = zeros(nSamples,1);
% 
% % generate a series of pulse time occurrences based on a Poisson process
% avStimRate = 100; % Hz, i.e. pulses per second
% dt = 1/avStimRate
% dS = round(fs * dt);
% lambda = dS; % lambda here pertains to the average stim rate, say 100 Hz
% r = poissrnd(lambda, 12000, 1);
% 
% % space out all the 1's according to the poisson process sample nums
% rCum = cumsum(r);
% rCum(rCum > nSamples) = [];
% stimEvents(rCum) = 1;
% 
% 
% % as a safeguard measure, control for all ISI's that would cause the pulse
% % waveform to continue on to the next one too
% 
% 
% % Convolve stim pulse waveform with event occurrences
% stimFile = conv(stimEvents, waveform); 
% lWav = length(waveform);
% stimFile(end-(lWav-2):end) = [];
% tFile = ((1:nSamples)-1) * (1/fs);
% figure; plot(tFile, stimFile)
