% code for generating stim file of variable IPI stim pulses for Lankarany
% collab




%% specify the waveform shape

srate = 24414.0625; % samples/sec; system sampling rate of TDT setup (e.g. 24414.0625)

dt = 1/srate;

zeropad = [0 0 0];
% Anode/Cathode length ratio:
aclr = 5;
catLen = 3;
waveformBase = [-aclr * ones(1,catLen), 0, ones(1,catLen*aclr)] / aclr;
% waveformBase = [-9 -9 -9 ones(1,27)];
% factor to multiply waveform by to get it to the desired current value
gain = 450; % microamps; 
waveform = [zeropad waveformBase zeropad] * gain;
t = 1e6 * [0:(length(waveform)-1)] * dt; % microseconds
figure; plot(t, waveform); xlabel('microsends'); ylabel('microamps')


%% Generate entire stim file with four components:
% 15 sec 130 Hz jittered
% 15 sec 130 Hz isochronous
% 15 sec 185 Hz jittered
% 15 sec 185 Hz isochronous


%   srate = 24414.0625; % samples/sec; system sampling rate of TDT setup (e.g. 24414.0625)
percStdms = 0.1; % scalar; the percentage around avISIms for 1 standard deviation (0.1),
% 10% of a milliseconds


% specify array of 4 segments to concatenate for stim conditions

totaltime = 15; % seconds; total time when stims will be on
%   dbsfreq = [130, 130, 185, 185]; % Hz
%   avISIms = [1e3*(1./(dbsfreq))]; % milliseconds; average inter-stimulus interval for the pulses

% stimIdx = [];
% 
% 
% 
% 
% for i = 1:4
%     i_stimIdx = genJitterStimtrain(totaltime(i), srate, avISIms(i), percStdms);
%     stimIdx = [stimIdx; i_stimIdx];
% 
% end


% create first stim file based on jittered 130Hz stim
dbsfreq = 130; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stimtrain(1).stimidx = genJitterStimtrain(totaltime, srate, avISIms, percStdms);
stimtrain(1).stimType = 'jittered';
stimtrain(1).srate = srate;
stimtrain(1).stimWav = waveform;
stimtrain(1).stimPulsepSec = dbsfreq;


% create first stim file based on isochronous 130Hz stim
dbsfreq = 130; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stimtrain(2).stimidx = genRegularStimtrain(totaltime, srate, avISIms, percStdms);
stimtrain(2).stimType = 'regular';
stimtrain(2).srate = srate;
stimtrain(2).stimWav = waveform;
stimtrain(2).stimPulsepSec = dbsfreq;

% create first stim file based on jittered 185Hz stim
dbsfreq = 185; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stimtrain(3).stimidx = genJitterStimtrain(totaltime, srate, avISIms, percStdms);
stimtrain(3).stimType = 'jittered';
stimtrain(3).srate = srate;
stimtrain(3).stimWav = waveform;
stimtrain(3).stimPulsepSec = dbsfreq;

% create first stim file based on isochronous 185Hz stim
dbsfreq = 185; % Hz
avISIms = 1e3*1/(dbsfreq); % milliseconds; average inter-stimulus interval for the pulses
stimtrain(4).stimidx = genRegularStimtrain(totaltime, srate, avISIms, percStdms);
stimtrain(4).stimType = 'regular';
stimtrain(4).srate = srate;
stimtrain(4).stimWav = waveform;
stimtrain(4).stimPulsepSec = dbsfreq;



stim = [stimtrain(1).stimidx; stimtrain(2).stimidx; ...
    stimtrain(3).stimidx; stimtrain(4).stimidx];

%% convolve the events with the stimulus waveform
stimFile = conv(stim, waveform); 
nSamples = length(stimFile);

lWav = length(waveform);
% stimFile(end-(lWav-2):end) = 0 ; % remove any partial waveform shapes
tFile = ((1:nSamples)-1) * (1/fs);
figure; plot(tFile, stimFile)

% save the results as a .mat file to be loaded into TDT Synapse software
save('stimtrain4cond_metadata.mat', 'stimtrain');
save('stimtrain4cond.mat', 'stimFile');






