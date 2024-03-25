% code for generating file stim in TDT Synapse for irregular DBS pulses
function [stim_OU] = genJitterStimtrain(totaltime, srate, avISIms, percStdms)
% calculates stim times in milliseconds
%
%
% INPUT
% total time: total time of stimtrain, seconds
% srate: sampling rate of system to deliver stim events, samples / second
% avISI: average interstimulus interval, milliseconds
% percStd: the percentage around avISIms for 1 standard deviation (0.1),
% 10% of a milliseconds
T = totaltime; % sec
% dt = 1e3*(1/srate); % msec
dt = (1/srate); % sec

Fs_DBS = 1/(avISIms/1e3); % Hz, pulses per sec

% All Milad's code below:
% T = 60; % sec
% dt = 0.1; % msec
% Fs_DBS = 100; % Hz
% mu_ISI = 1e3*1/Fs_DBS; % (10 msec eq. = 100 Hz), msec
mu_ISI = 1/Fs_DBS; % (0.01 sec eq. = 100 Hz), sec
std_ISI = percStdms * mu_ISI; % sec

stim = zeros(floor(T/dt * 1e3),1);

threshold_min_ISI = 2/1000;% sec the timing of two consecutive pulses cannot be less than 2 msec
num_pulse = T*Fs_DBS;
% rand_ISI = randn(num_pulse,1);
% ISI = mu_ISI + std_ISI*rand_ISI; % msec
% ISI(ISI<=threshold_min_ISI) = threshold_min_ISI;
% stim(floor(cumsum(ISI/dt))) = 1;
% stim(floor(T/dt * 1e3):end) = 0; 
% stim = stim(1:floor(T/dt * 1e3));

% figure; plot(dt:T*1e3/dt,max(ISI)*stim)
%hold on, plot((1:length(ISI))/dt/1e3,ISI,'b')

% figure; plot((1:length(ISI))/dt/1e3,ISI,'k')


%% USE OU Process
% tau = 1; % msec
% ISI_OU = ISI; 
% % inp_sig = zeros(size(ISI)); 
% inp_sig = zeros(size(ISI)+1); 
% 
% ETA = inp_sig;
% N_tau = sqrt(2/tau);
% k=1;
% for t = dt:dt:T/dt
%     inp_sig(k+1) = rand_ISI(k);
%     Einf = tau*N_tau*inp_sig(k)/sqrt(dt);
%     ETA(k + 1) = Einf + (ETA(k) - Einf)*exp(-dt/tau); 
%     k = k+1;
% end
% ISI_OU = mu_ISI + std_ISI*ETA;
% ISI_OU(ISI_OU<=threshold_min_ISI) = threshold_min_ISI;
% 
% stim_OU = zeros(floor(T/dt * 1e3),1);
% stim_OU(floor(cumsum(ISI_OU/dt))) = 1;
% stim_OU(floor(T/dt * 1e3):end) = 0; 
% stim_OU = stim_OU(1:floor(T/dt * 1e3));

% figure; plot(dt:T*1e3/dt,max(ISI_OU)*stim_OU)
%hold on, plot((1:length(ISI))/dt/1e3,ISI,'b')

% figure; plot((1:length(ISI_OU))/dt/1e3,ISI_OU,'k')


%% My final version of using the Orhstein-Uhlenbeck process in Matlab


th = 1;
mu = 0;
sig = 1;
% dt = 1e-2;
% dt = mu_ISI;
t = 0:mu_ISI:T;             % Time vector
% t = 0:dt:T;             % Time vector
x0 = 0;                 % Set initial condition
rng(2);                 % Set random seed
W = zeros(1,length(t)); % Allocate integrated W vector
for i = 1:length(t)-1
    W(i+1) = W(i)+sqrt(exp(2*th*t(i+1))-exp(2*th*t(i)))*randn;
end
ex = exp(-th*t);
x = x0*ex+mu*(1-ex)+sig*ex.*W/sqrt(2*th);



ISI_OU = mu_ISI + std_ISI*x;
% figure; histogram(ISI_OU); figure; plot(t,ISI_OU);
any(isinf(ISI_OU))
any(isnan(ISI_OU))
ISI_OU(ISI_OU<=threshold_min_ISI) = threshold_min_ISI;

% clean up from any inf and nan values
ISI_OU(isnan(ISI_OU)) = [];
ISI_OU(isinf(ISI_OU)) = [];


stim_OU = zeros(floor(T*srate),1);
stim_OU(floor(cumsum(ISI_OU*srate))) = 1; % assign 1's for the stim time indices
% stim_OU(floor(T/dt * 1e3):end) = 0; 
% stim_OU = stim_OU(1:floor(T/dt * 1e3));



end