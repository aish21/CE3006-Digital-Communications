clear all; 
close all; 
clc;

%---------------Define Declaration---------------%
%Num of bits = 1024 (given)
nBits = 1024;

%assuming signal power to be one as given               
Signal_Power = 1;  

%SNR_dB = 10log(Signal_Power/Noise_Power)
%Create matrix to store SNR
SNR_dB = 0:0.5:20; 
%might need to change (Consider different SNR values from 0 dB to 50 dB (in multiples of 5 dB))

%SNR = S/N = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10)); %converting to the normal SNR 

%Holder value for plotting later
plot_signal = rand(1,nBits);
plot_noise = rand(1,nBits);
plot_receive = rand(1,nBits);

%MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
%TO PLOT SIGNAL,NOISE and RECEIVE
plot_SNR_dB = 15;

%Counter for number of run to calculate BER for each SNR value
Total_Run = 20;
%---------------END OF DEFINE--------------------%


%----------------MAIN ROUTINE--------------------%

%Iterate through different SNR value
for i = 1 : length(SNR)
	Sum_Error = 0;
	for j = 1 : Total_Run
        % for each snr run 20 times
        
        %-----------TRANSMITTER------------%
        
		%Generate random binary digits(0 or 1) of length 1024, INPUT SIGNAL
		Data = round(rand(1,nBits)); 
        
		%Convert binary digit to -1 or +1: (-/+ 1 coding given)
        %if Data is 0 -> signal = 2*-0.5 = -1
		%else if Data is 1 -> signal = 2*0.5 = +1
		Signal = 2 .* (Data - 0.5);
        
		%Get Noise Power from SNR (SNR = SP/NP)
		Noise_Power = Signal_Power ./SNR(i);
        
        %Get Noise Variance from Noise Power
        %Change the noise variance with respect to SNR (signal to noise ratio) value
	%sigma square = N0/2
        Avg_Noise_Power = Noise_Power/2;
        
        %The randn() function is for normal distribution, generated with equal number of noise samples 
        %(In Matlab, a.*randn(1000,1) + b is a Gaussian Dist, mean b = 0 sd a = root(variance N0)) 
        Noise = sqrt(Avg_Noise_Power) .*randn(1,nBits);   
        
		%Receiver side Signal (noise is added to the signal)
		Receive = Signal+Noise;                      

        
        %-----RECEIVER------%
        
        %Initialize threshold
		Threshold = 0; %(the transmitted data is +1 and -1, and 0 is the mid value)
        
        %Initialize Error for this run
		Error = 0;
        Output = zeros(1,nBits); %array of zeros (one row, 1024 cols)
        %Fix the threshold value as 0
        %If received signal >= threshold value, threshold = 1
        %If received signal < threshold value, threshold = 0
		for k= 1 : nBits
            if (Receive(k)>= Threshold)
                Output(k) = 1;
            end
            if (Receive(k) < Threshold)
                Output(k) = 0;
            end
			if (Receive(k)>= Threshold) && Data(k)==0||(Receive(k)<Threshold && Data(k)==1) %i think can just compare output and data, instead of using threshold again???
				
                Error = Error+1;
			end
        end
        
        
        %---------------ERROR STATISTICS CALCULATION--------------------% 
		%BER = TotalError divide by number of bits
		Error = Error ./nBits;  
        
		%Accumulate the error of each run within the SNR	
		Sum_Error = Error + Sum_Error;
        
        %---------------Choose the value for plotting------------------%
        if (plot_SNR_dB == SNR_dB(i))
            plot_signal = Signal;
            plot_noise = Noise;
            plot_receive = Receive;
            plot_output = Output;
        end
        
    end
    
    %Average Error for that particular SNR
	Error_Rate(i) = Sum_Error / Total_Run;
end

%Predict BER using 
Pred_BER=(1/2)*erfc(sqrt(SNR)); 


%--------------------------------PLOTTING--------------------------------%
figure("position", [10,100,1400,800]) 

%BER against SNR -- main plot
figure(1);
semilogy (SNR_dB,Error_Rate,'r*');
xlabel('Normalized SNR')
ylabel('Probability Error');
title('BER against SNR');
hold on
semilogy (SNR_dB,Pred_BER,'m');
legend('Simulation','Prediction');
axis([0 20 10^(-5) 1]);
hold off

%data generation
figure(2)
subplot(221);plot(plot_signal);title('Data Generated');
subplot(224);plot(plot_noise);title('Noise Generated');
subplot(222);plot(plot_receive);title('Received Data');
subplot(223);plot(plot_output);title('Output');



