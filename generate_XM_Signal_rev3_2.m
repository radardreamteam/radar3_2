%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate the XM signal (QPSK modulation) 
%
%   All Eq. references can be found in the following paper:
%       'Use of XM radio satellite signal as a source of opportunity for
%       passive coherent location' by L.P. Gill, et. al.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all 

generatePlot = 1;
%% Define Constants
j = 1i;
propSpeed 	= 299792458; % m/s

%%Simulation parameters

directSigPower = -122.94+30;
echoSigPower   = -193.93+30;
refGain        = 20;
dirPathAttenuation                 = -50;%dB
survAntGain                        = 40;

%% System Parameters
SystemParameters;

%% Create input signal
% sigNumber   = 400;
pilotOn     = 1; % insert pilot module
numFrames   = 1;
% sigLength   = sampsPerCycle*cyclesPerSymbol*sigNumber;
sigLength   = DVBLength(cyclesPerSymbol,sampsPerCycle,numFrames,pilotOn);
dirPath     = zeros(1, sigLength + samp_offset);
indirPath   = zeros(1, sigLength + samp_offset);
tempEcho = zeros(1, sigLength + samp_offset);
tempDir  =  zeros(1, sigLength + samp_offset );
sig = zeros(1, sampsPerCycle*cyclesPerSymbol*400+samp_offset );
taxis = (0:(length(dirPath)-1))/samplingFreq;

%% Time domain signal
%______________________________QPSK random_________________________________
% dirPath(1:sigLength)               = ZgenrateQPSK_Signal(cyclesPerSymbol,sampsPerCycle, fcSignal,sigNumber);
%______________________________DVB signal_________________________________
pilotOn = 1; % insert pilot module
dirPath(1:sigLength)               = DVBgenerate(cyclesPerSymbol,sampsPerCycle, fcSignal,numFrames,pilotOn);

indirPath(samp_offset+1:end)       = dirPath(1:sigLength)*exp(j*phaseOffset).*exp(j*2*pi*FShift*taxis(1:sigLength));
tempEcho(samp_offset+1:end)        = indirPath(samp_offset+1:end); 

%The direct signal from the Tx should arrive at both Rx at the same time,
%there is no shift.
tempDir                            = dirPath;
%make the reference signal
refSignal                          = 10^((directSigPower+refGain)/10)*tempDir;
%Create surveillance channel
% dirPathAttenuation                 = -50;%dB
% survAntGain                        = 40;
echoSigPower                       = echoSigPower+survAntGain;
echoSignal                         = 10^(echoSigPower/10)*tempEcho;
directSignal                       = 10^(dirPathAttenuation/10)*10^(directSigPower/10)*tempDir;
survChannel                        = echoSignal +directSignal ;

%add noise
noiseSig                           = 10^(noisePower_dBm/10)*(randn(1,length(tempDir))+j*randn(1,length(tempDir)));
survNoisyChannel = survChannel +noiseSig;

noiseSig                           = 10^(noisePower_dBm/10)*(randn(1,length(tempDir))+j*randn(1,length(tempDir)));
NoisyrefSignal = refSignal+noiseSig;
freqVector = -300:2:300;

%-------------------------------------------------------------------------
%nlms
filterOrder = 32;
nlms = dsp.LMSFilter(filterOrder,'Method','Normalized LMS','StepSizeSource','Input port');

[y,err,weights] = nlms(NoisyrefSignal',survNoisyChannel',0.001);


%-------------------------------------------------------------------------
[rdmap, ranges, freqs] = rangedopplerfft(survNoisyChannel',samplingFreq , 2*timeDelay*propSpeed , freqVector, NoisyrefSignal');
[rdmap_nlms, ranges_nlms, freqs_nlms] = rangedopplerfft(err,samplingFreq , 2*timeDelay*propSpeed , freqVector, NoisyrefSignal');
[X,Y] = meshgrid(ranges, freqs);

f1 = figure('Name','QPSK random'); 
contourf(X,Y,rdmap');
savefig('.\myPlot\QPSK_random.fig');

f2 = figure('Name','QPSK random after NLMS'); 
contourf(X,Y,rdmap_nlms');
saveas(f2,'.\myPlot\QPSK_random_after_NLMS.fig')

% Frequency domain signal
dirPath_fft = fftshift(abs(fft(dirPath))/length(dirPath));
indirPath_fft = fftshift(abs(fft(indirPath))/length(indirPath));

%% Process Data


%% Plots
% axes
taxis_sig = (0:(length(dirPath)-1))/samplingFreq;
faxis_sig = linspace(-samplingFreq/2,samplingFreq/2,length(dirPath));

% Time domain of signals
f3 = figure; 
hold on; plot(taxis_sig*1e6, real(dirPath));
hold on; plot(taxis_sig*1e6, real(indirPath)); axis tight; grid
xlabel('Time (\mus)')
ylabel('Amplitude (V)')
title('Time Domain XM Signal')
legend('Direct Path','Indirect Path')
savefig(f3,'.\myPlot\Signals.fig')

% Freq domain of signals
f4 = figure(); 
indFreqSide = find((faxis_sig >= 0),1);

hold on; plot(faxis_sig(indFreqSide:end)*1e-6, 20*log10(dirPath_fft(indFreqSide:end)));

axis([0 10 -200 5]); grid
xlabel('Frequency (MHz)')
ylabel('Amplitude (dBm)')
title('XM signal Spectrum')
%legend('XM signal','Direct Path','Indirect Path')
saveas(f4,'.\myPlot\Signal_spectrum','fig');


%% DOA estimation using MUSIC algorithm
% Build a 15*15 rectangular array, and use MUSIC algorithm to seperate these two signals 
array = phased.URA('Size',[15 15],'ElementSpacing',[lamda/2 lamda/2]);
array.Element.FrequencyRange = [500.0e6 5000.0e6];
doa1 = [25;10];
doa2 = [30;-20];
survChannelArray  = collectPlaneWave(array,[echoSignal',directSignal'],[doa1,doa2],samplingFreq);
noiseArray        = 10^(noisePower_dBm/10)*(randn(length(tempDir),225)+j*randn(length(tempDir),225));
survChannelArray  = survChannelArray + 0*noiseArray ;
%% NLMS
%We must assume we know the doa of the direct path signal
direct_distribute = collectPlaneWave(array,directSignal',doa2,samplingFreq);
noiseArray        = 10^(noisePower_dBm/10)*(randn(length(tempDir),225)+j*randn(length(tempDir),225));
direct_distribute = direct_distribute+noiseArray;
NewsurvChannelArray = zeros(size(survChannelArray,1),size(survChannelArray,2));
parfor i = 1:size(survChannelArray,2)
    nlms.reset();
    [y,err,weights] = nlms(direct_distribute(:,i),survChannelArray(:,i),0.001);
    NewsurvChannelArray(:,i) = err;
end


%% range-doppler map 
%compensate the phase shift to each element
reference_distribute = collectPlaneWave(array,NoisyrefSignal',doa2,samplingFreq);
multi_channel_rdmap = zeros(size(rdmap,1),size(rdmap,2),size(reference_distribute,2));
parfor i = size(reference_distribute,2)
    [rdmap, ranges, freqs] = rangedopplerfft(NewsurvChannelArray(:,i),samplingFreq , 2*timeDelay*propSpeed , freqVector, reference_distribute(:,i));
    multi_channel_rdmap(:,:,i) = rdmap;
end
rdmap_compansated = sum(multi_channel_rdmap,3);


f5 = figure('Name','Phased array QPSK random after NLMS and compensation'); 
contourf(X,Y,rdmap_compansated')
saveas(f5,'.\myPlot\Array_NLMS_compensated','fig')

analogBeam = sum(NewsurvChannelArray,2);
[rdmap, ranges, freqs] = rangedopplerfft(analogBeam,samplingFreq , 2*timeDelay*propSpeed , freqVector, NoisyrefSignal');
f6 = figure('Name','Phased array QPSK random after NLMS without compensation'); 
contourf(X,Y,rdmap')
saveas(f6,'.\myPlot\Array_NLMS_not_compensated','fig')

% hrdresp = phased.RangeDopplerResponse(...
%    'RangeMethod','FFT',...
%    'PropagationSpeed',propSpeed,...
%    'SampleRate',samplingFreq,...
%    'DechirpInput',false,...
%    'SweepSlope',RangeDopplerEx_Dechirp_SweepSlope);
% [resp,rng_grid,dop_grid] = response(RangeDopplerEx_MF_X, ...
%     RangeDopplerEx_MF_Coeff);

%% angle estimation
estimator = phased.MUSICEstimator2D('SensorArray',array,...
    'OperatingFrequency',samplingFreq,...
    'NumSignalsSource','Property',...
    'DOAOutputPort',true,'NumSignals',2,...
    'AzimuthScanAngles',-30:.5:30,...
    'ElevationScanAngles',-30:.5:30);
[~,doas] = estimator(NewsurvChannelArray)
f7 = figure(); 
plotSpectrum(estimator);
title('MUSIC algorithm');
saveas(f7,'.\myPlot\Array_NLMS_not_compensated','fig');



%% Bistatic parameter estimation
 [tx2tgEst,rx2tgEst] = paraExtraction(Ddiff, norm(tx2rx), theta1);



