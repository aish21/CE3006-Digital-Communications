clear all; 
clear workspace;

% Number of encoded bits
nBits = 1792;
% Given carrier frequency
carrierFrequency = 10000; 
% Sampled at 16 times the frequency
samplingFrequency = 32 * carrierFrequency;
% Given data rate
dataRate = 1000;

% Generate the random signal
Signal = randi([0 1],1,1024);
% Doing hamming encoding
hammingSignal = encode(Signal,7,4,'hamming/binary');
% Max timestamp
totalTime = nBits/dataRate;
% Calculate sampling period
samplingPeriod = 1/samplingFrequency;
% Calculate number of data points
nSteps = totalTime/samplingPeriod;
% Defining the time points
t = 0:samplingPeriod:totalTime;

% Define the 2 carrier functions
carrier1 = cos(2*pi*carrierFrequency*t);
carrier0 = cos(pi*carrierFrequency*t);

% Thresholds for the 6th order butterworth filter
[b,a] = butter(6,0.2);

% Find the transmitted signal
transmittedSignal = zeros([1 nSteps]);
count = 0;
% Calculating the number of samples per bit of signal
timePerBit = 1/dataRate;
noOfSamplesPerBit = timePerBit / samplingPeriod;
% Creating the sampled signal
for k = 1: nSteps
    transmittedSignal(k) = hammingSignal(ceil(k*dataRate/samplingFrequency));
end
transmittedSignal(k + 1) = transmittedSignal(k);

% Find the modulated signal
FSKmodulated1 = carrier1 .* (transmittedSignal == 1);
FSKmodulated2 = carrier0 .* (transmittedSignal == -1);
modulatedSig = FSKmodulated1 + FSKmodulated2;
% Calculating the signal power
sigPower = rms(modulatedSig)^2;

index = 0;
% Defining the SNR range
SNRvalues = -50:5:50;
meanBitError = zeros([1 length(SNRvalues)]);
% Iterating through the different SNR values
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
        % decoding the hamming signal
        finalResultDecoded = decode(result, 7, 4, 'hamming/binary');
        % Calculating error
        bitErrorRate(sample) = mean(finalResultDecoded ~= Signal);
    end

    % Assigning to the mean error rate
    meanBitError(index) = mean(bitErrorRate);

    % Variable Plots
    if (SNR == 5)
        plotFSK = modulatedSig;
        plotFSKTransmit = transmittedFSK;
        plotLPfilterFSK = BFSK_demod;
        plotfinalResultFSK = result;
    end
end

% Plotting the mean BER vs SNR
figure(1)
hold on
semilogy(SNRvalues, meanBitError,'k-*');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('FSK');

hold off

% Plotting the different stages of FSK modulation and demodulation
figure(2);
titles = {'Data Generated: ', 'Step - Modulation (FSK): ', 'Signal Recieved: ','Step - Demodulation (FSK): ','Decoded: '}; 
subplot(5,1,1);plot(Signal);title('Data Generated: ');
subplot(5,1,2);plot(plotFSK,'k');title('Step - Modulation (FSK): ');
subplot(5,1,3);plot(plotFSKTransmit, 'k');title('Signal Recieved: ');
subplot(5,1,4);plot(plotLPfilterFSK, 'k');title('Step - Demodulation (FSK): ');
subplot(5,1,5);plot(plotfinalResultFSK);title('Decoded: ');
