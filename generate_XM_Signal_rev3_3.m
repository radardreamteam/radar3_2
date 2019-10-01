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
sigNumber   = 400; 
pilotOn     = 1; % insert pilot module
numFrames   = 1;
threshold = 1e-15;
debug = false;


for intIdx = 1:length(sigNumber)
    
    sigLength   = sampsPerCycle*cyclesPerSymbol*sigNumber(intIdx);
    %sigLength   = DVBLength(cyclesPerSymbol,sampsPerCycle,numFrames,pilotOn);
    dirPath     = zeros(1, sigLength + samp_offset);
    indirPath   = zeros(1, sigLength + samp_offset);
    tempEcho = zeros(1, sigLength + samp_offset);
    tempDir  =  zeros(1, sigLength + samp_offset );
    sig = zeros(1, sampsPerCycle*cyclesPerSymbol*400+samp_offset );
    taxis = (0:(length(dirPath)-1))/samplingFreq;

    %% Time domain signal
    %______________________________QPSK random_________________________________
    dirPath(1:sigLength)               = ZgenrateQPSK_Signal(cyclesPerSymbol,sampsPerCycle, fcSignal,sigNumber(intIdx));
    %dirPath(1:sigLength)               = DVBgenerate(cyclesPerSymbol,sampsPerCycle, fcSignal,numFrames,pilotOn);
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


    figure(1);
    subplot(1, length(sigNumber), intIdx);
    contourf(X,Y,rdmap');
    [centroidsscaled1,height1,prob_det1]=compute_centroids(X,Y,rdmap',threshold,debug);
    hold on;
    plot(centroidsscaled1(1),centroidsscaled1(2),'b*');
    hold off;
    xlabel('Range (m)');
    ylabel('Doppler (m/s)');
    title(['QPSK random (IntTime = ' num2str(sigNumber(intIdx)) ')']);

    figure(2);
    subplot(1, length(sigNumber), intIdx);
    contourf(X,Y,rdmap_nlms');
    [centroidsscaled2,height2,prob_det2]=compute_centroids(X,Y,rdmap_nlms',threshold,debug);
    hold on;
    plot(centroidsscaled2(1),centroidsscaled2(2),'b*');
    hold off;
    xlabel('Range (m)');
    ylabel('Doppler (m/s)');
    title(['QPSK random after NLMS (IntTime = ' num2str(sigNumber(intIdx)) ')']);


    % Frequency domain signal
    dirPath_fft = fftshift(abs(fft(dirPath))/length(dirPath));
    indirPath_fft = fftshift(abs(fft(indirPath))/length(indirPath));

    %% Process Data


    %% Plots
    % axes
    taxis_sig = (0:(length(dirPath)-1))/samplingFreq;
    faxis_sig = linspace(-samplingFreq/2,samplingFreq/2,length(dirPath));

    % Time domain of signals
    figure(3);
    subplot(1, length(sigNumber), intIdx);
    hold on; plot(taxis_sig*1e6, real(dirPath));
    hold on; plot(taxis_sig*1e6, real(indirPath)); axis tight; grid
    xlabel('Time (\mus)')
    ylabel('Amplitude (V)')
    title(['Time Domain XM Signal (IntTime = ' num2str(sigNumber(intIdx)) ')'])
    legend('Direct Path','Indirect Path')


    % Freq domain of signals
    figure(4);
    subplot(1, length(sigNumber), intIdx);
    indFreqSide = find((faxis_sig >= 0),1);

    hold on; plot(faxis_sig(indFreqSide:end)*1e-6, 20*log10(dirPath_fft(indFreqSide:end)));

    axis([0 10 -200 5]); grid
    xlabel('Frequency (MHz)')
    ylabel('Amplitude (dBm)')
    title(['XM signal Spectrum (IntTime = ' num2str(sigNumber(intIdx)) ')'])
    %legend('XM signal','Direct Path','Indirect Path')



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


    figure(5);
    subplot(1, length(sigNumber), intIdx);
    contourf(X,Y,rdmap_compansated')
    [centroidsscaled3,height3,prob_det3]=compute_centroids(X,Y,rdmap_compansated',threshold,debug);
    hold on;
    plot(centroidsscaled3(1),centroidsscaled3(2),'b*');
    hold off;
    xlabel('Range (m)');
    ylabel('Doppler (m/s)');
    title(['Phased array QPSK random after NLMS with compensation (IntTime = ' num2str(sigNumber(intIdx)) ')']);


    analogBeam = sum(NewsurvChannelArray,2);
    [rdmap, ranges, freqs] = rangedopplerfft(analogBeam,samplingFreq , 2*timeDelay*propSpeed , freqVector, NoisyrefSignal');
    figure(6);
    subplot(1, length(sigNumber), intIdx);
    contourf(X,Y,rdmap')
    [centroidsscaled4,height4,prob_det4]=compute_centroids(X,Y,rdmap',threshold,debug);
    hold on;
    plot(centroidsscaled4(1),centroidsscaled4(2),'b*');
    hold off;
    xlabel('Range (m)');
    ylabel('Doppler (m/s)');
    title(['Phased array QPSK random after NLMS without compensation (IntTime = ' num2str(sigNumber(intIdx)) ')']);


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
    figure(7); 
    subplot(1, length(sigNumber), intIdx);
    plotSpectrum(estimator);
    title(['MUSIC algorithm (IntTime = ' num2str(sigNumber(intIdx)) ')']);


end

%% Save Figures - Comment out sections whether using Windows or Mac/Linux

% For Windows
% saveas(figure(1),'.\myPlot\QPSK_random.fig', 'fig');
% saveas(figure(2),'.\myPlot\QPSK_random_after_NLMS.fig', 'fig')
% saveas(figure(3),'.\myPlot\Signals.fig', 'fig')
% saveas(figure(4),'.\myPlot\Signal_spectrum.fig','fig');
% saveas(figure(5),'.\myPlot\Array_NLMS_compensated','fig')
% saveas(figure(6),'.\myPlot\Array_NLMS_not_compensated','fig')
% saveas(figure(7),'.\myPlot\MUSIC_algorithm','fig');

% % For Mac/Linux
% saveas(figure(1),'./myPlot/QPSK_random.fig', 'fig');
% saveas(figure(2),'./myPlot/QPSK_random_after_NLMS.fig', 'fig')
% saveas(figure(3),'./myPlot/Signals.fig', 'fig')
% saveas(figure(4),'./myPlot/Signal_spectrum.fig','fig');
% saveas(figure(5),'./myPlot/Array_NLMS_compensated','fig')
% saveas(figure(6),'./myPlot/Array_NLMS_not_compensated','fig')
% saveas(figure(7),'./myPlot/MUSIC_algorithm','fig');

%% Bistatic parameter estimation
 [tx2tgEst,rx2tgEst] = paraExtraction(Ddiff, norm(tx2rx), theta1);



