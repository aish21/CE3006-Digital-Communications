close all;
clear workspace;
clc;

% Modulation Technique 1: Binary Phase Shift Keying 

nBits = 1024;                               %Num of bits = 1024 (given)
carrierFrequency = 10000;                   %Let the carrier frequency be 10KHz
samplingFreq = carrierFrequency * 16;       %signal is 16 times oversampled 
dataRate = 1000;                            %considering the baseband data rate as 1 kbps.
samplingPeriod = samplingFreq / dataRate;   
[lowB, lowA] = butter(6,0.2);               %6th order filter with cut-off frequency 0.2 (low pass)

%Define Amplitude, time
amp = 1;                                    
t = 0: 1/samplingFreq : nBits/dataRate;

SNR_dB = -50:5:50;                         %Initialize SNR range and steps
SNR = (10.^(SNR_dB/10));                   %Converting SNR DB
errorRateBPSK = zeros(1,length(SNR));      %Intialize error variable
signalLength = samplingFreq*nBits/dataRate + 1; %calculate singal length after sampling

carrier = amp .* cos(2*pi*carrierFrequency*t); %define carrier signal

numRuns = 10;                                  %define number times to run 

%Change to inspect waveforms for different SNRs
SelectedSNR_DB = 5;

for i = 1 : length(SNR)
    avgError = 0;
    for j = 1 : numRuns
        data = round(rand(1,nBits));
        dataStream = zeros(1, signalLength);
        
        for k = 1: signalLength - 1
            dataStream(k) = data(ceil(k*dataRate/samplingFreq));
        end
        dataStream(signalLength) = dataStream(signalLength - 1);
    
        bpskSourceSignal = 2.* dataStream - 1;                      %Make the signal 1 and -1
        bpskSignal = carrier .* bpskSourceSignal; 
           
        %calculate expected noise for transmission  
        bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
        bpskNoisePower = bpskSignalPower ./ SNR(i);
        bpskNoise = sqrt(bpskNoisePower/2) .* randn(1,signalLength);
        
        %adding noise to the generated signal
        bpskTransmit = bpskSignal + bpskNoise;
        
        % Use coherent detection - to demodulate transmitted signal
        bpskDemodulated = bpskTransmit .* (2 .* carrier);
        
        % Pass this through the LP filter
        bpskFiltered = filtfilt(lowB, lowA,  bpskDemodulated);
        
        %Sample the signal
        bpskSample = sample(bpskFiltered, samplingPeriod, nBits);
        finalResultBPSK = decision_device(bpskSample,nBits, 0);
    
        %Calculate Error 
        errorBPSK = 0;
        for k = 1: nBits - 1
            if(finalResultBPSK(k) ~= data(k))
                errorBPSK = errorBPSK + 1;
            end
        end 
        avgError = avgError + errorBPSK/nBits;
        
        % Variable Plots
        if (SelectedSNR_DB == SNR_dB(i))
            plotSignal = data;
            plotBPSK = bpskSignal;
            plotBPSKTransmit = bpskTransmit;
            plotLPfilterBPSK = bpskFiltered;
            plotfinalResultBPSK = finalResultBPSK;
        end
    end
    errorRateBPSK(i) = avgError/numRuns;
end

%Plotting the BER vs SNR graph
figure(1);
plot1 = semilogy (SNR_dB,errorRateBPSK,'k-*');
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
title("Bit error probability for BPSK modulation")
hold off

%plot waveforms at different stages 
figure(2);
subplot(511);plot(plotSignal);title('Data Generated: ');
subplot(512);plot(plotBPSK,'k');title('Step - Modulation (BPSK): ');
subplot(513);plot(plotBPSKTransmit, 'k');title('Signal Recieved: ');
subplot(514);plot(plotLPfilterBPSK, 'k');title('Step - Demodulation (BPSK): ');
subplot(515);plot(plotfinalResultBPSK);title('Decoded: ');

% Sampling and Decision Device Simulation
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
