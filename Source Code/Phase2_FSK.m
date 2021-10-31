clear all; 
clear workspace;

% Given number of bits
nBits = 1024;
% Given carrier frequency
carrierFrequency = 10000; 
% Sampled at 16 times the frequency
samplingFrequency = 32 * carrierFrequency;
% Given data rate
dataRate = 1000;

% Generate the random signal
Signal = randi([0 1],1,nBits);
% Max timestamp
totalTime = nBits/dataRate;
% Calculate the time instances
samplingPeriod = 1/samplingFrequency;
nSteps = totalTime/samplingPeriod;
t = 0:samplingPeriod:totalTime;

% Define the 2 carrier functions
carrier1 = cos(2*pi*carrierFrequency*t);
carrier0 = cos(pi*carrierFrequency*t);

% Thresholds for the 6th order butterworth filter
[b,a] = butter(6,0.2);

% Find the transmitted signal
transmittedSignal = zeros([1 nSteps]);
count = 0;
timePerBit = 1/dataRate;
noOfSamplesPerBit = timePerBit / samplingPeriod;
for bit = Signal
    for i = 1:noOfSamplesPerBit
        count = count+1;
        transmittedSignal(count) = bit;
    end
end
transmittedSignal(count + 1) = transmittedSignal(count);

% Find the modulated signal
FSKmodulated1 = carrier1 .* (transmittedSignal == 1);
FSKmodulated2 = carrier0 .* (transmittedSignal == -1);
modulatedSig = FSKmodulated1 + FSKmodulated2;
sigPower = rms(modulatedSig)^2;

index = 0;
SNRvalues = -50:5:50;
meanBitError = zeros([1 length(SNRvalues)]);
for SNR = SNRvalues
    index = index+1;
    bitErrorRate = zeros([1 nBits]);
    % Calculate noise power from SNR
    noisePower = sigPower/(10^(SNR/10));

    % Run the experiment for each sample 20 times
    for sample = 1:20
        % Generate the noise
        Noise = randn(1,length(modulatedSig));
        Noise = sqrt(noisePower) .* Noise;
        
        % Find the transmitted signal after noise
        transmittedFSK = modulatedSig+Noise;

        %Coherent demodulation of FSK
        BFSK_demod1 = transmittedFSK*2.* carrier1;
        BFSK1_filter = filtfilt(b,a,BFSK_demod1);
        BFSK_demod0 = transmittedFSK*2 .* carrier0;
        BFSK0_filter = filtfilt(b,a,BFSK_demod0);
        BFSK_demod = BFSK_demod1 -BFSK_demod0 ;

        count = 0;
        result = zeros([1 nBits]);
        % sampling
        for i = 20:noOfSamplesPerBit:length(BFSK_demod)
            count = count+1;
            result(count) = BFSK_demod(i);
        end

        % Applying threshold logic
        for i = 1:length(result)
            if result(i) >= 0.5
                result(i) = 1;
            else
                result(i) = 0;
            end
        end

        % Calculating error
        bitErrorRate(sample) = mean(result ~= Signal);
    end
    meanBitError(index) = mean(bitErrorRate);

    % Variable Plots
    if (SNR == 5)
        plotFSK = modulatedSig;
        plotFSKTransmit = transmittedFSK;
        plotLPfilterFSK = BFSK_demod;
        plotfinalResultFSK = result;
    end
end
% 

figure(1)
%semilogy (SNRvalues, theoreticalError,'r', 'linewidth', 1.5);
hold on
semilogy(SNRvalues, meanBitError,'k-*');
hold off

figure(2);
subplot(5,1,1);title('Data Generated: ');plot(Signal);
subplot(5,1,2);title('Step - Modulation (FSK): ');plot(plotFSK,'k');
subplot(5,1,3);title('Signal Recieved: ');plot(plotFSKTransmit, 'k')
subplot(5,1,4);title('Step - Demodulation (FSK): ');plot(plotLPfilterFSK, 'k');
subplot(5,1,5);title('Decoded: ');plot(plotfinalResultFSK);
