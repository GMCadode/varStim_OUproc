% script for generating a charge-balanced stim pulse waveform with
% cathodic-anodic symmetric shape. cathodic part is ~120 microseconds long,
% same for anodic, with total pulse time of ~240 microseconds. 



fs = 24414.0625; % sampling frequency that I expect to use in TDT
dt = 1/fs;

zeropad = [0 0 0];
waveformBase = [-1 -1 -1 1 1 1];
gain = 100/1e6; % factor to multiply waveform by to get it to the desired current value

waveform = [zeropad waveformBase zeropad] * gain;
t = [0:(length(waveform)-1)] * dt;
figure; plot(t, waveform)
