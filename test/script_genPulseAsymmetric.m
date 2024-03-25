% script for generating a charge-balanced stim pulse waveform, with
% cathodic-assymetric shape. cathodic part is ~120 microseconds long. Ratio
% of cathodic/anodic amplitude is 9/1 (90% to 10%). Total pulse duration is
% ~1.4 milliseconds.



fs = 24414.0625; % sampling frequency that I expect to use in TDT
dt = 1/fs;

zeropad = [0 0 0];
waveformBase = [-9 -9 -9 ones(1,27)];
% factor to multiply waveform by to get it to the desired current value
gain = (100/1e6) / 9; % in this case, 100 microAmps
waveform = [zeropad waveformBase zeropad] * gain;
t = [0:(length(waveform)-1)] * dt;
figure; plot(t, waveform)
