%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jordan Hayes
%Project 1: Automated Event Detection 
%EGR434: Bioelectric Potentials 
%Dr. Rhodes 
%10/28/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Import data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ECG_dataTable = readtable("EGR434_Project1_RestingandExerciseEKG_jordan.txt");
%ECG_dataTable = readtable("Project1_data2.txt");
t = ECG_dataTable(:,1); 
t = table2array(t);

ECG = ECG_dataTable(:,4);
ECG = table2array(ECG);
ECG(isnan(ECG)) = 0;

ECG_resting = ECG(1:2001,1);
ECG_postExercise = ECG(2002:4001,1); 


n = length(ECG);
N = 30;

%zero mean 
%ECG = ECG - mean(ECG); 
ECG_resting = ECG_resting - mean(ECG_resting); 
ECG_postExercise = ECG_postExercise - mean(ECG_postExercise); 

ECG = zeros(4001,1); 
ECG(1:2001, 1) = ECG_resting(:,1); 
ECG(2002:4001, 1) = ECG_postExercise(:,1);

ECG = ECG - mean(ECG); 
%Normalize to 1
ECG = ECG ./ max(abs(ECG));

Ts = t(2,1) - t(1,1);                           %numerical diferrence in indices = Ts 
fs = 1 / Ts;                                    %1/Ts = fs 

figure(1)
plot(t, ECG); 
grid on;
xlabel('t (s)'); 
ylabel('ECG (uV)'); 
title('ECG vs Time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 1, Low Pass Filter
% LPF (1-z^-6)^2/(1-z^-1)^2
%given in Pans Tompkins study
% Eq 1, where z^-x denotes a weight of one in index x
% and the constant is the value in the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];                         %numermator
a=[1 -2 1];                                             %denominator 
 
 
LP_filter = filter(b,a,[1, zeros(1,length(b))]);        % transfer function of LPF

%filter signal(s) with Low Pass
ECG_step1 = conv(ECG ,LP_filter);

%filter has a shift of 6
ECG_step1 = ECG_step1(6+[1: n]);                         %cancel delay

%normalize again 
ECG_step1 = ECG_step1 ./ max(abs(ECG_step1)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting for test block 
% figure(2)
% plot([0:length(ECG_step1) - 1]/fs, ECG_step1); 
% xlim([0, 20])
% grid on;
% xlabel('t (s)'); 
% ylabel('ECG (uV)'); 
% title('LP Filtered ECG vs Time'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 2, High Pass Filter
% HPF = Pass - (Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
%given in Pans Tompkins study
% Eq 4, where z^-x denotes a weight of one in index x
% and the constant is the value in the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High pass filter 

b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];
 
HP_filter = filter(b,a,[1, zeros(1,length(b))]);                %impulse response iof HPF

%filter signal with HP
ECG_step2 = conv(ECG_step1, HP_filter);

%filter processing delay of 16
ECG_step2 = ECG_step2(16+[1: n]);

%normalize
ECG_step2 = ECG_step2 ./ max(abs(ECG_step2));

% figure(3)
% plot([0:length(ECG_step2) - 1]/fs, ECG_step2); 
% xlim([0, 20])
% grid on;
% xlabel('t (s)'); 
% ylabel('ECG (uV)'); 
% title('HP Filtered ECG vs Time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 3, Derivative Filter
% HPF = Pass - (Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
%given in Pans Tompkins study
% Eq 7, where z^-x denotes a weight of one in index x
% and the constant is the value in the index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%derivative filter = transfer function impulse response
h = [-1 -2 0 2 1]/8;                                                    %impulse response given in Pans Tomp Study

%apply derivative filter
ECG_step3 = conv(ECG_step2, h); 

%shift delay
ECG_step3 = ECG_step3(2+[1:n]);

%normalize
ECG_step3 = ECG_step3 ./ max(abs(ECG_step3));
 
% figure(4)
% plot([0:length(ECG_step3) - 1]/fs, ECG_step3); 
% xlim([0, 20])
% grid on;
% xlabel('t (s)'); 
% ylabel('ECG (uV)'); 
% title('Derivative Filtered ECG vs Time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 5. Squaring
% y(nT) = [x(nT)]^2 which means index times sampling period
%indicating a point in time. same thing as x(t)

ECG_step4 = ECG_step3 .^ 2; 

%normalize
ECG_step4 = ECG_step4 ./ max(abs(ECG_step4));

% figure(5)
% plot([0:length(ECG_step4) - 1]/fs, ECG_step4); 
% xlim([0, 20])
% grid on;
% xlabel('t (s)'); 
% ylabel('ECG (uV)'); 
% title('Squared ECG vs Time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 6. Moving Window Integration 
% given in Pans Tompkins study, 
% Eq. 11
% The purpose of moving-window integration is to 
% obtain waveform feature information in addition 
% to the slope of the R wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For our sample rate of 200 samples/s, 
% the window is 30 samples wide (150 ms)

h = ones(1 ,N+1)/N;                 %n = 30, indices 1:31

%apply filter
ECG_step5 = conv(ECG_step4, h); 

%cancel delay of 15 samples 
ECG_step5 = ECG_step5 (15+[1: n]);

%normalize
ECG_step5 = ECG_step5 ./ max(abs(ECG_step5)); 


%Plotting for test block 
% figure(6)
% plot([0:length(ECG_step5) - 1]/fs, ECG_step5); 
% xlim([0, 20])
% grid on;
% xlabel('t (s)'); 
% ylabel('ECG (uV)'); 
% title('Integration Window ECG vs Time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%detect QRS 

%set initial values to apply
pk1 = max(ECG_step5);
thresh = mean(ECG_step5);

%find peaks
[pks, pk_indices] = findpeaks(ECG_step5, 'MinPeakHeight', thresh*pk1, 'MinPeakDistance', fs / 5); 
 

ispk = (ECG_step5 > pk1 * thresh)';                             %find where V > thresh exists, sets binary vector

rising_edge = find(diff([0 ispk]) == 1);                        %will find the rising edges find(ispk - 0) == 1 where it hits threshold from 0 to 1
falling_edge = find(diff([ispk 0]) == -1);                      %will find the falling edges, or when the signal falls back below threshold
                                                  
pk_count = length(pks);

for i = 1:pk_count 


     [R_mag(i) R_index(i)] = max(ECG(rising_edge(i):falling_edge(i)));               %will assume the peak in time domain us R peak
     R_index(i) = R_index(i)-1+rising_edge(i);                                       % add offset
      
     [Q_mag(i) Q_index(i)] = min(ECG(rising_edge(i):R_index(i)));                    %under the rising edge is the window for QRS 
     Q_index(i) = Q_index(i)-1+rising_edge(i);                                       %q is to left of R, so min(rising, Rpk)
     
     [S_mag(i) S_index(i)] = min(ECG(rising_edge(i):falling_edge(i)));               %same logic but to the right, assuming lowest point
     S_index(i) = S_index(i)-1+rising_edge(i); 
     
     QRS_width(i) = (S_index(i) - Q_index(i)) * Ts; 
     
end

isRestingPk = find(pk_indices < 2002)';                     %find peaks during resting  
resting_pks = length(isRestingPk); 

isPostExercisePk = find(pk_indices > 2001)';                %find peaks during post exercise   
postExercise_pks = length(isPostExercisePk);                %length of the array is number of peaks


resting_heartRate = resting_pks * 6;                        %number of beats per minute
postExercise_heartRate = postExercise_pks * 6;              %number of events / minutes (10 sec = 1/6 min)

total_avg_heartRate = (resting_pks + postExercise_pks) * 3; %beats/min of total 20 second signal 

fprintf("The average resting heartrate is %0.2f\n", resting_heartRate);
fprintf("The average post exercise heartrate is %0.2f,\n", postExercise_heartRate);
fprintf("The overall average heartrate is %0.2f\n", total_avg_heartRate);

figure(7)
title('ECG Signal with R points');
%plot(t, ECG, t(R_index), R_mag, 'r^'); 
plot(t, ECG, t(R_index), R_mag, 'r^', t(S_index), S_mag, '*', t(Q_index), Q_mag, 'o');
legend('ECG','R','S','Q');
xlabel('t (s)'); 
ylabel('ECG (uV)'); 
title('ECG vs Time with QRS Points Idetified'); 

figure(8)
histogram(QRS_width);
title("QRS Width Distribution"); 
ylabel("Frequency")
xlabel("QRS Width (s)"); 




