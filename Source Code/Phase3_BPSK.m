close all;
clear workspace;
clc;

% Modulation Technique 1: Binary Phase Shift Keying 

nBits = 1024;                               %Num of bits = 1024 (given)
codeword_len = 7;                           %length of codeword
msgword_len = 4;                            %length of message word
encBits = nBits*codeword_len/msgword_len;   %no. of bits in encoded signal

%initialisations for CBC encoding
genpoly = cyclpoly(codeword_len,msgword_len); %n=7, k=4
parmat = cyclgen(codeword_len,genpoly);
trt = syndtable(parmat);

carrierFrequency = 10000;                   %Let the carrier frequency be 10KHz
samplingFreq = carrierFrequency * 16;       %signal is 16 times oversampled 
dataRate = 1000;                            %considering the baseband data rate as 1 kbps.
samplingPeriod = samplingFreq / dataRate;   
[lowB, lowA] = butter(6,0.2);               %6th order filter with cut-off frequency 0.2 (low pass)

%Define Amplitude, time
amp = 1;                                    
t = 0: 1/samplingFreq : nBits/dataRate;
tEnc = 0: 1/samplingFreq : encBits/dataRate;
SNR_dB = -50:5:50;                         %Initialize SNR range and steps
SNR = (10.^(SNR_dB/10));                   %Converting SNR DB
errorRateBPSK = zeros(1,length(SNR));      %Intialize error variable
errorRateBPSKHam = zeros(1,length(SNR));      %Intialize error variable for encoding HAMMING
errorRateBPSKCBC = zeros(1,length(SNR));      %Intialize error variable for encoding CBC

signalLength = samplingFreq*nBits/dataRate + 1; %calculate singal length after sampling
signalLengthenc = samplingFreq*encBits/dataRate + 1; %calculate singal length after sampling for encoding

carrier = amp .* cos(2*pi*carrierFrequency*t); %define carrier signal
carrierEnc = amp .* cos(2*pi*carrierFrequency*tEnc); %carrier signal for encoding

numRuns = 10;                                  %define number times to run 

%Change to inspect waveforms for different SNRs
SelectedSNR_DB = 5;

for i = 1 : length(SNR)
    avgError = 0;
    avgErrorHam = 0;
    avgErrorCBC = 0;

    for j = 1 : numRuns
        data = round(rand(1,nBits));
        dataStream = zeros(1, signalLength);
        
        %hamming method initialisation
        hammingSignal = encode(data,codeword_len,msgword_len,'hamming/binary');
        dataStreamHam = zeros(1, signalLengthenc);
        
        cbcSignal = encode(data,codeword_len,msgword_len,'cyclic/binary',genpoly);
        dataStreamCBC = zeros(1, signalLengthenc);
        
        for k = 1: signalLength - 1
            dataStream(k) = data(ceil(k*dataRate/samplingFreq));
        end
        dataStream(signalLength) = dataStream(signalLength - 1);
        bpskSourceSignal = 2.* dataStream - 1;                      %Converts zeros to -1 and 1 remains 1
        bpskSignal = carrier .* bpskSourceSignal;                   %bpskSignal is the modulated signal
        
        for k = 1: signalLengthenc - 1
            dataStreamHam(k) = hammingSignal(ceil(k*dataRate/samplingFreq));
            dataStreamCBC(k) = cbcSignal(ceil(k*dataRate/samplingFreq));
        end
        dataStreamHam(signalLengthenc) = dataStreamHam(signalLengthenc - 1);
        bpskSourceSignalHam = 2.* dataStreamHam - 1;                      %Converts zeros to -1 and 1 remains 1
        bpskSignalHam = carrierEnc .* bpskSourceSignalHam;                %bpskSignalHam is the modulated signal for HAMMING

        dataStreamCBC(signalLengthenc) = dataStreamCBC(signalLengthenc - 1);
        bpskSourceSignalCBC = 2.* dataStreamCBC - 1;                      %Converts zeros to -1 and 1 remains 1
        bpskSignalCBC = carrierEnc .* bpskSourceSignalCBC;                %bpskSignalCBC is the modulated signal for CBC
           
        %calculate expected noise for transmission  
        bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
        bpskNoisePower = bpskSignalPower ./ SNR(i);
        bpskNoise = sqrt(bpskNoisePower/2) .* randn(1,signalLength);
        
        %calculate expected noise for transmission - HAMMING
        bpskSignalPowerHam = (norm(bpskSignalHam)^2)/signalLengthenc;
        bpskNoisePowerHam = bpskSignalPowerHam ./ SNR(i);
        bpskNoiseHam = sqrt(bpskNoisePowerHam/2) .* randn(1,signalLengthenc);

        %calculate expected noise for transmission - CBC
        bpskSignalPowerCBC = (norm(bpskSignalCBC)^2)/signalLengthenc;
        bpskNoisePowerCBC = bpskSignalPowerCBC ./ SNR(i);
        bpskNoiseCBC = sqrt(bpskNoisePowerCBC/2) .* randn(1,signalLengthenc);
        
        %adding noise to the generated signal
        bpskTransmit = bpskSignal + bpskNoise;
        
        %adding noise to the generated signal - HAMMING
        bpskTransmitHam = bpskSignalHam + bpskNoiseHam;

        %adding noise to the generated signal - CBC
        bpskTransmitCBC = bpskSignalCBC + bpskNoiseCBC;
        
        % Use coherent detection - to demodulate transmitted signal
        bpskDemodulated = bpskTransmit .* (2 .* carrier);
        
        % Use coherent detection - to demodulate transmitted signal - HAMMING
        bpskDemodulatedHam = bpskTransmitHam .* (2 .* carrierEnc);
        
        % Use coherent detection - to demodulate transmitted signal - CBC
        bpskDemodulatedCBC = bpskTransmitCBC .* (2 .* carrierEnc);

        % Pass this through the LP filter
        bpskFiltered = filtfilt(lowB, lowA,  bpskDemodulated);
        
        % Pass this through the LP filter - HAMMING
        bpskFilteredHam = filtfilt(lowB, lowA,  bpskDemodulatedHam);

        % Pass this through the LP filter - CBC
        bpskFilteredCBC = filtfilt(lowB, lowA,  bpskDemodulatedCBC);
        
        %Sample the signal
        bpskSample = sample(bpskFiltered, samplingPeriod, nBits);
        finalResultBPSK = decision_device(bpskSample,nBits, 0);
        
        %Sample the signal - HAMMING
        bpskSampleHam = sample(bpskFilteredHam, samplingPeriod, encBits);
        finalResultBPSKHam = decision_device(bpskSampleHam,encBits, 0);
        finalResultBPSKDecodedHam = decode(finalResultBPSKHam, codeword_len, msgword_len, 'hamming/binary');
    
        %Sample the signal - CBC
        bpskSampleCBC = sample(bpskFilteredCBC, samplingPeriod, encBits);
        finalResultBPSKCBC = decision_device(bpskSampleCBC,encBits, 0);
        finalResultBPSKDecodedCBC = decode(finalResultBPSKCBC,codeword_len,msgword_len,'cyclic/binary',genpoly,trt);
    
        %Calculate Error 
        errorBPSK = 0;
        errorBPSKHam = 0;
        errorBPSKCBC=0;
        
        for k = 1: nBits - 1
            %Calculate Error normal
            if(finalResultBPSK(k) ~= data(k))
                errorBPSK = errorBPSK + 1;
            end
            %Calculating Error HAMMING
            if(finalResultBPSKDecodedHam(k) ~= data(k))
                errorBPSKHam = errorBPSKHam + 1;
            end
            %Calculating Error CBC
            if(finalResultBPSKDecodedCBC(k) ~= data(k))
                errorBPSKCBC = errorBPSKCBC + 1;
            end
        end 
        avgError = avgError + errorBPSK/nBits;
        avgErrorHam = avgErrorHam + errorBPSKHam/nBits;
        avgErrorCBC = avgErrorCBC + errorBPSKCBC/nBits;
        
        % Variable Plots
        if (SelectedSNR_DB == SNR_dB(i))
            plotSignal = data;
            plotBPSK = bpskSignal;
            plotBPSKTransmit = bpskTransmit;
            plotLPfilterBPSK = bpskFiltered;
            plotfinalResultBPSK = finalResultBPSK;
            
            %HAMMING
            plotSignalHam = hammingSignal;
            plotBPSKHam = bpskSignalHam;
            plotBPSKTransmitHam = bpskTransmitHam;
            plotLPfilterBPSKHam = bpskFilteredHam;
            plotfinalResultBPSKHam = finalResultBPSKDecodedHam;

            %CBC
            plotSignalCBC = cbcSignal;
            plotBPSKCBC = bpskSignalCBC;
            plotBPSKTransmitCBC = bpskTransmitCBC;
            plotLPfilterBPSKCBC = bpskFilteredCBC;
            plotfinalResultBPSKCBC = finalResultBPSKDecodedCBC;
        end
    end
    errorRateBPSK(i) = avgError/numRuns;
    errorRateBPSKHam(i) = avgErrorHam/numRuns;
    errorRateBPSKCBC(i) = avgErrorCBC/numRuns;
end

%Plotting the BER vs SNR graph
figure(1);
semilogy (SNR_dB,errorRateBPSK,'r');
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
title("BPSK - Encoding Comparison")
hold on

semilogy (SNR_dB,errorRateBPSKHam,'b');
hold on

semilogy (SNR_dB,errorRateBPSKCBC,'g');
leg = legend('BPSK','BPSK-Hamming','BPSK-CBC');
set(leg,'location','southwest')
hold off

%plot waveforms at different stages 
figure(2);
subplot(511);plot(plotSignal);title('Data Generated: ');
subplot(512);plot(plotBPSK,'k');title('Step - Modulation (BPSK): ');
subplot(513);plot(plotBPSKTransmit, 'k');title('Signal Recieved: ');
subplot(514);plot(plotLPfilterBPSK, 'k');title('Step - Demodulation (BPSK): ');
subplot(515);plot(plotfinalResultBPSK);title('Decoded: ');

%plot waveforms at different stages HAMMING
figure(3);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalHam);title('Encoded Data - Hamming: ');
subplot(613);plot(plotBPSKHam,'k');title('Step - Modulation (BPSK): ');
subplot(614);plot(plotBPSKTransmitHam, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterBPSKHam, 'k');title('Step - Demodulation (BPSK): ');
subplot(616);plot(plotfinalResultBPSKHam);title('Decoded: ');

%plot waveforms at different stages CBC
figure(4);
subplot(611);plot(plotSignal);title('Data Generated: ');
subplot(612);plot(plotSignalCBC);title('Encoded Data - CBC: ');
subplot(613);plot(plotBPSKCBC,'k');title('Step - Modulation (BPSK): ');
subplot(614);plot(plotBPSKTransmitCBC, 'k');title('Signal Recieved: ');
subplot(615);plot(plotLPfilterBPSKCBC, 'k');title('Step - Demodulation (BPSK): ');
subplot(616);plot(plotfinalResultBPSKCBC);title('Decoded: ');

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
