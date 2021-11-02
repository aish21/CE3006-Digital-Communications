clear all; 
close all;
clear workspace;

% Number of encoded bits
nBits = 1024;                               %Num of bits = 1024 (given)
codeword_len = 7;                           %length of codeword
msgword_len = 4;                            %length of message word
encBits = nBits*codeword_len/msgword_len;   %no. of bits in encoded signal

%initialisations for CBC encoding
genpoly = cyclpoly(codeword_len,msgword_len); 
parmat = cyclgen(codeword_len,genpoly);
trt = syndtable(parmat);

% Given carrier frequency
carrierFrequency = 10000; 
% Sampled at 16 times the frequency
samplingFrequency = 16 * carrierFrequency;
% Given data rate
dataRate = 1000;

% Generate the random signal
Signal = randi([0 1],1,1024);
% Doing hamming encoding
hammingSignal = encode(Signal,codeword_len,msgword_len,'hamming/binary');
% Doing cyclic block code encoding
cbcSignal = encode(Signal,codeword_len,msgword_len,'cyclic/binary',genpoly);
% Max timestamp
totalTime = nBits/dataRate;
totalEncTime = encBits/dataRate;
% Calculate sampling period
samplingPeriod = 1/samplingFrequency;
% Calculate number of data points
nSteps = totalTime/samplingPeriod;
nStepsEnc = totalEncTime/samplingPeriod;
% Defining the time points
t = 0:samplingPeriod:totalTime;
tEnc = 0:samplingPeriod:totalEncTime;

% Define the 2 carrier functions
carrier1 = cos(2*pi*carrierFrequency*t);
carrier0 = cos(pi*carrierFrequency*t);

% Define the 2 carrier functions for encoded signals
EncCarrier1 = cos(2*pi*carrierFrequency*tEnc);
EncCarrier0 = cos(pi*carrierFrequency*tEnc);

% Thresholds for the 6th order butterworth filter
[b,a] = butter(6,0.2);

% Find the transmitted signal
transmittedSignal = zeros([1 nSteps]);
transmittedHamSignal = zeros([1 nStepsEnc]);
transmittedCBCSignal = zeros([1 nStepsEnc]);
count = 0;

% Calculating the number of samples per bit of signal
timePerBit = 1/dataRate;
noOfSamplesPerBit = timePerBit / samplingPeriod;

% Creating the sampled signal
for k = 1: nSteps
    transmittedSignal(k) = Signal(ceil(k*dataRate/samplingFrequency));
end
transmittedSignal(k + 1) = transmittedSignal(k);

for k = 1: nStepsEnc
    transmittedHamSignal(k) = hammingSignal(ceil(k*dataRate/samplingFrequency));
    transmittedCBCSignal(k) = cbcSignal(ceil(k*dataRate/samplingFrequency));
end
transmittedHamSignal(k + 1) = transmittedHamSignal(k);
transmittedCBCSignal(k + 1) = transmittedCBCSignal(k);

% Find the modulated signal
FSKmodulated1 = carrier1 .* (transmittedSignal == 1);
FSKmodulated2 = carrier0 .* (transmittedSignal == -1);
modulatedSig = FSKmodulated1 + FSKmodulated2;

%Hamming
FSKHammodulated1 = EncCarrier1 .* (transmittedHamSignal == 1);
FSKHammodulated2 = EncCarrier0 .* (transmittedHamSignal == -1);
modulatedHamSig = FSKHammodulated1 + FSKHammodulated2;

%CBC
FSKCBCmodulated1 = EncCarrier1 .* (transmittedCBCSignal == 1);
FSKCBCmodulated2 = EncCarrier0 .* (transmittedCBCSignal == -1);
modulatedCBCSig = FSKCBCmodulated1 + FSKCBCmodulated2;

% Calculating the signal power
sigPower = rms(modulatedSig)^2;
HamsigPower = rms(modulatedHamSig)^2;
CBCsigPower = rms(modulatedCBCSig)^2;

index = 0;
% Defining the SNR range
SNRvalues = -50:5:50;
meanBitError = zeros([1 length(SNRvalues)]);
meanBitErrorHam = zeros([1 length(SNRvalues)]);
meanBitErrorCBC = zeros([1 length(SNRvalues)]);
% Iterating through the different SNR values
for SNR = SNRvalues
    index = index+1;
    bitErrorRate = zeros([1 nBits]);
    bitErrorRateHam = zeros([1 nBits]);
    bitErrorRateCBC = zeros([1 nBits]);
    % Calculate noise power from SNR
    noisePower = sigPower/(10^(SNR/10));
    noisePowerHam = HamsigPower/(10^(SNR/10));
    noisePowerCBC = CBCsigPower/(10^(SNR/10));

    % Run the experiment for each sample 20 times
    for sample = 1:20
        % Generate the noise
        Noise = randn(1,length(modulatedSig));
        Noise = sqrt(noisePower) .* Noise;

        % Generate the noise for Hamming encoded Signal
        HammingNoise = randn(1,length(modulatedHamSig));
        HammingNoise = sqrt(noisePowerHam) .* HammingNoise;


        % Generate the noise for CBC encoded Signal
        CBCNoise = randn(1,length(modulatedCBCSig));
        CBCNoise = sqrt(noisePowerCBC) .* CBCNoise;
        
        % Find the transmitted signal after noise
        transmittedFSK = modulatedSig+Noise;
        transmittedFSKHam = modulatedHamSig+HammingNoise;
        transmittedFSKCBC = modulatedCBCSig+CBCNoise;

        %Coherent demodulation of FSK
        BFSK_demod1 = transmittedFSK*2.* carrier1;
        BFSK1_filter = filtfilt(b,a,BFSK_demod1);
        BFSK_demod0 = transmittedFSK*2 .* carrier0;
        BFSK0_filter = filtfilt(b,a,BFSK_demod0);
        BFSK_demod = BFSK1_filter -BFSK0_filter ;

        %Coherent demodulation of FSK - HAMMING
        BFSK_demod1_Ham = transmittedFSKHam*2.* EncCarrier1;
        BFSK1_filter_Ham = filtfilt(b,a,BFSK_demod1_Ham);
        BFSK_demod0_Ham = transmittedFSKHam*2 .* EncCarrier0;
        BFSK0_filter_ham = filtfilt(b,a,BFSK_demod0_Ham);
        BFSK_demod_ham = BFSK1_filter_Ham - BFSK0_filter_ham ;

        %Coherent demodulation of FSK - CBC
        BFSK_demod1_CBC = transmittedFSKCBC*2.* EncCarrier1;
        BFSK1_filter_CBC = filtfilt(b,a,BFSK_demod1_CBC);
        BFSK_demod0_CBC = transmittedFSKCBC*2 .* EncCarrier0;
        BFSK0_filter_CBC = filtfilt(b,a,BFSK_demod0_CBC);
        BFSK_demod_CBC = BFSK1_filter_CBC -BFSK0_filter_CBC ;

        % sampling
        count = 0;
        result = zeros([1 nBits]);
        for i = 20:noOfSamplesPerBit:length(BFSK_demod)
            count = count+1;
            result(count) = BFSK_demod(i);
        end

        count = 0;
        resultHam = zeros([1 encBits]);
        % sampling HAMMING
        for i = 20:noOfSamplesPerBit:length(BFSK_demod_ham)
            count = count+1;
            resultHam(count) = BFSK_demod_ham(i);
        end

        %Sampling CBC
        resultCBC = zeros([1 encBits]);
        count = 0;
        for i = 20:noOfSamplesPerBit:length(BFSK_demod_CBC)
            count = count+1;
            resultCBC(count) = BFSK_demod_CBC(i);
        end

        % Applying threshold logic
        for i = 1:length(result)
            if result(i) >= 0.5
                result(i) = 1;
            else
                result(i) = 0;
            end
        end

        % Applying threshold logic - HAMMING
        for i = 1:length(resultHam)
            if resultHam(i) >= 0.5
                resultHam(i) = 1;
            else
                resultHam(i) = 0;
            end
        end

        % Applying threshold logic - CBC
        for i = 1:length(resultCBC)
            if resultCBC(i) >= 0.5
                resultCBC(i) = 1;
            else
                resultCBC(i) = 0;
            end
        end
        
        % decoding the hamming signal
        finalResultDecodedHam = decode(resultHam, codeword_len,msgword_len, 'hamming/binary');
        % decoding the cbc signal
        finalResultDecodedCBC = decode(resultCBC, codeword_len,msgword_len,'cyclic/binary',genpoly,trt);
        
        % Calculating error
        bitErrorRate(sample) = mean(result ~= Signal);
        bitErrorRateHam(sample) = mean(finalResultDecodedHam ~= Signal);
        bitErrorRateCBC(sample) = mean(finalResultDecodedCBC ~= Signal);
    end

    % Assigning to the mean error rate
    meanBitError(index) = mean(bitErrorRate);
    meanBitErrorHam(index) = mean(bitErrorRateHam);
    meanBitErrorCBC(index) = mean(bitErrorRateCBC);

    % Variable Plots
    if (SNR == 5)
        plotSignal = Signal;
        plotFSK = modulatedSig;
        plotFSKTransmit = transmittedFSK;
        plotLPfilterFSK = BFSK_demod;
        plotfinalResultFSK = result;
        
        %HAMMING
        plotSignalHam = hammingSignal;
        plotFSKHam = modulatedHamSig;
        plotFSKTransmitHam = transmittedFSKHam;
        plotLPfilterFSKHam = BFSK_demod_ham;
        plotfinalResultFSKHam = finalResultDecodedHam;

        %CBC
        plotSignalCBC = cbcSignal;
        plotFSKCBC = modulatedCBCSig;
        plotFSKTransmitCBC = transmittedFSKCBC;
        plotLPfilterFSKCBC = BFSK_demod_CBC;
        plotfinalResultFSKCBC = finalResultDecodedCBC;
    end
end

% Plotting the mean BER vs SNR
figure(1)
semilogy(SNRvalues, meanBitError,'r');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('FSK - Encoding Comparison');
hold on

semilogy(SNRvalues, meanBitErrorHam,'b');
hold on

semilogy(SNRvalues, meanBitErrorCBC,'g');
leg = legend('FSK','FSK-Hamming','FSK-CBC');
set(leg,'location','southwest')
hold off

% Plotting the different stages of FSK modulation and demodulation
figure(2);
titles = {'Data Generated: ', 'Step - Modulation (FSK): ', 'Signal Recieved: ','Step - Demodulation (FSK): ','Decoded: '}; 
subplot(5,1,1);plot(plotSignal);title('Data Generated: ');
subplot(5,1,2);plot(plotFSK,'k');title('Step - Modulation (FSK): ');
subplot(5,1,3);plot(plotFSKTransmit, 'k');title('Signal Recieved: ');
subplot(5,1,4);plot(plotLPfilterFSK, 'k');title('Step - Demodulation (FSK): ');
subplot(5,1,5);plot(plotfinalResultFSK);title('Decoded: ');

%plot waveforms at different stages HAMMING
figure(3);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalHam);title('Encoded Data - Hamming: ');
subplot(613);plot(plotFSKHam,'k');title('Step - Modulation (FSK): ');
subplot(614);plot(plotFSKTransmitHam, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterFSKHam, 'k');title('Step - Demodulation (FSK): ');
subplot(616);plot(plotfinalResultFSKHam);title('Decoded: ');

%plot waveforms at different stages CBC
figure(4);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalCBC);title('Encoded Data - CBC: ');
subplot(613);plot(plotFSKCBC,'k');title('Step - Modulation (FSK): ');
subplot(614);plot(plotFSKTransmitCBC, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterFSKCBC, 'k');title('Step - Demodulation (FSK): ');
subplot(616);plot(plotfinalResultFSKCBC);title('Decoded: ');
