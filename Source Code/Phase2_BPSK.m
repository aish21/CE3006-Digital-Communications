close all;
clc;

% Modulation Technique 1: Binary Phase Shift Keying 

%{ 
    - Num of bits = 1024 (given)
    - Let the carrier frequency be 10KHz
    - Assume one cycle contains 16 samples
    - Baseband data rate = 1Kbps
    - Assume a 6th order filter with cut-off frequency 0.2 (low pass)
    - Define Amplitude, time
    - Initialize SNR, set value to change plots
    - Initialize error variable
    - Calculate signal length
    - Set number of runs
    - Functions to simulate decision device + sampling
    - Plotting
%}

nBits = 1024; %Num of bits = 1024 (given)
carrierFrequency = 10000; %Let the carrier frequency be 10KHz
samplingFreq = carrierFrequency * 16; % signal is oversampled 16 time
dataRate = 1000;
samplingPeriod = samplingFreq / dataRate;
[lowB, lowA] = butter(6,0.2);
amp = 1;
t = 0: 1/samplingFreq : nBits/dataRate;
SNR_dB = -50:5:50;
SNR = (10.^(SNR_dB/10));
errorRateBPSK = zeros(1,length(SNR));
carrier = amp .* cos(2*pi*carrierFrequency*t);
signalLength = samplingFreq*nBits/dataRate + 1;
numRuns = 10;

modifySNR_dB = 5;

for i = 1 : length(SNR)
    avgError = 0;
    for j = 1 : numRuns
        data = round(rand(1,nBits));
        dataStream = zeros(1, signalLength);
        
        for k = 1: signalLength - 1
            dataStream(k) = data(ceil(k*dataRate/samplingFreq));
        end
        dataStream(signalLength) = dataStream(signalLength - 1);
    
        bpskSourceSignal = 2.* dataStream - 1; %*2 bc of negative 1 and 1 bc ook is 0 to 1 so don't need *2
        bpskSignal = carrier .* bpskSourceSignal; % this too
           
        bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
        bpskNoisePower = bpskSignalPower ./ SNR(i);
        bpskNoise = sqrt(bpskNoisePower/2) .* randn(1,signalLength);
    
        bpskTransmit = bpskSignal + bpskNoise;
    
        bpskDemodulated = bpskTransmit .* (2 .* carrier);
        bpskFiltered = filtfilt(lowB, lowA,  bpskDemodulated);
    
        bpskSample = sample(bpskFiltered, samplingPeriod, nBits);
        finalResultBPSK = decision_device(bpskSample,nBits, 0);
    
        errorBPSK = 0;
       
        for k = 1: nBits - 1
            %fprintf("%d -> %d\n", data(k), finalResultBPSK(k));
            if(finalResultBPSK(k) ~= data(k))
                errorBPSK = errorBPSK + 1;
            end
        end 
        avgError = avgError + errorBPSK/nBits;
        
        % Variable Plots
        if (modifySNR_dB == SNR_dB(i))
            plotSignal = data;
            plotBPSK = bpskSignal;
            plotBPSKTransmit = bpskTransmit;
            plotLPfilterBPSK = bpskFiltered;
            plotfinalResultBPSK = finalResultBPSK;
        end
    end
    errorRateBPSK(i) = avgError/numRuns;
end

figure(1);
plot1 = semilogy (SNR_dB,errorRateBPSK,'k-*');
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
title("BPSK")
hold off

figure(2);
subplot(511);title('Data Generated: ');plot(plotSignal);
subplot(512);title('Step - Modulation (OOK): ');plot(plotBPSK,'k');
subplot(513);title('Signal Recieved: ');plot(plotBPSKTransmit, 'k')
subplot(514);title('Step - Demodulation (OOK): ');plot(plotLPfilterBPSK, 'k');
subplot(515);title('Decoded: ');plot(plotfinalResultBPSK);

function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end

function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end
