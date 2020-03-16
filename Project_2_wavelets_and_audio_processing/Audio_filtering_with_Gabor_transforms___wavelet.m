%%%%%%%%Audio Processing with wavelets%%%%%%%%%
%%%%%%%%%         Josef Henthorn     %%%%%%%%%
%%%%%%%%%          Winter 2020       %%%%%%%%% 
%%%%%%%%%   University of Washington %%%%%%%%%
%
% Objective: 
%   1. Filter audio signal data 
%   2. Use gabor temporal filtering
%   3. Use wavelets for temporal filtering
%   4. Build spectrograms and other useful graphs
%
%
% SUMMARY: The following code is a demonstration of ... 
%  Filtering a signal in the temporal domain and conducting tme-frequency
%  analysis.
% 
clear; close all; clc;
tic

load handel
v = y';

% p8 = audioplayer(v,Fs);
% playblocking(p8);
 

L = (1/Fs)*length(v); % Seconds: 1/(sample rate) * # of samples
n = length(v);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
b = (1:length(v))/Fs;
tau = 100; % Width of Gabor/Gaussian filter
samp_step = 0.01;
gabor = exp(-tau*(t-4.5).^2); % Gabor/Gaussian filter
gabor_2 = exp(-tau*(t-(4.5-(L*samp_step))).^2); % Gabor/Gaussian filter
gabor_3 = exp(-tau*(t-(4.5+(L*samp_step))).^2); % Gabor/Gaussian filter

vg = gabor.*v; % Gabor filtered signal
vgt = abs(fftshift(fft(vg))); % Transformed gabor signal

k = (2*pi/(L))*[0:(n/2) -n/2:-1]; ks = fftshift(k); % Frequency domain
vt = abs(fftshift(fft(v)));



% figure(1) = figure('Position', [200, 100, 1250, 1250]);
% subplot(3,5,[1 2])
% plot(b,v,'b'), hold on
% plot(b,gabor,'r','Linewidth',3)
% %plot(b,gabor_2,'r','Linewidth',3)
% %plot(b,gabor_3,'r','Linewidth',3)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of handle audio and Gabor filter (\tau = 0.1)');
% set(gca,'Fontsize',[12])
% 
% subplot(3,5,[6 7])
% plot((1:length(v))/Fs,vg)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title({'Gabor filtered signal of handle audio'});
% set(gca,'Fontsize',[12])
% 
% subplot(3,5,[11 12])
% plot(ks,vgt/max(vgt))
% xlabel('\omega [Hz]');
% ylabel('Amplitude');
% title({'Frequency spectrum of Gabor filtered handle audio'});
% set(gca,'Fontsize',[12])
%saveas(gcf,'3 graphs: tau = 100, samp step: 0.08925 sec shift.png')
%print('3 graphs: tau = 10','-djpeg')

%%
% GAUSSIAN FUNCTION
%gaussian = exp( -((t - 4.5).^2 ));
% MEXICAN HAT FUNCTION
%mexican_hat = (1-(t-4.5).^2).*exp((-(t-4.5).^2)/2); % Mexican Hat filter
% UNIT STEP FUNCTION
%unitstep = zeros(size(t)); 
%unitstep(6>=t & t>=4) = 1;

% figure(2) = figure('Position', [200, 100, 1250, 1250]);
% subplot(1,3,1)
% plot(t,gaussian,'k','Linewidth',3)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% grid on
% title('Gaussian wavelet filter');
% set(gca,'Fontsize',[15])
% 
% subplot(1,3,2)
% plot(t,mexican_hat,'k','Linewidth',3)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% grid on
% title('Mexican Hat wavelet filter');
% set(gca,'Fontsize',[15])
% 
% subplot(1,3,3)
% plot(t,unitstep,'k','Linewidth',3)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% grid on
% title('Shannon wavelet filter');
% set(gca,'Fontsize',[15])
% 
% print('Gabor wavelet filters','-djpeg')
%%
samp_step = 0.01;
tau = 100; % Width of Gabor/Gaussian filter

% Gaussian filtering for loop
tslide = 0:samp_step:L;

vgt_spec = [];
signal_mexihat_spec = [];
unitstep = zeros(size(t));
signal_shannon_spec = [];

for j = 1:length(tslide)
  
    gab = exp(-tau*(t-tslide(j)).^2); % Gabor/Gaussian filter
    vg = gab.*v; % Gabor filtered signal
    vgt = abs(fftshift(fft(vg))); % Transformed gabor signal
    vgt_spec = [vgt_spec; vgt/max(vgt)];
    
    mexican_hat = (1-(t-tslide(j)).^2).*exp((-tau*(t-tslide(j)).^2)); % Mexican Hat filter
    filter_sig_mh = mexican_hat.*v;
    filter_sig_mx_fft = abs(fftshift(fft(filter_sig_mh)));
    signal_mexihat_spec = [signal_mexihat_spec; filter_sig_mx_fft/max(filter_sig_mx_fft)];
    
    %unitstep((t-tslide(j) & t+tslide(j))) = 1;
    shannon = heaviside(t-tslide(j)) - heaviside(t-(1+tslide(j)));
    filter_sig_shannon = shannon.*v;
    filter_sig_shannon_fft = abs(fftshift(fft(filter_sig_shannon)));
    signal_shannon_spec = [signal_shannon_spec; filter_sig_shannon_fft/max(filter_sig_shannon_fft)];
    

      
%     filter_sig_s = unitstep.*v;
%     filter_sig_s_fft = abs(fftshift(fft(filter_sig_s)));
%     signal_shannon_spec = [signal_shannon_spec; filter_sig_s_fft/max(filter_sig_s_fft)];
%     figure(3)
%     subplot(3,1,1)
%     plot(b,v,'b'), hold on
%     plot(b,shannon,'r','Linewidth',3), hold off
%     xlabel('Time [sec]');
%     ylabel('Amplitude');
%     title('Signal of Interest with Gaussian filter');    
    
    
    
%     figure(3)
%     subplot(3,1,1)
%     plot(b,v,'b'), hold on
%     plot(b,gab,'r','Linewidth',3), hold off
%     xlabel('Time [sec]');
%     ylabel('Amplitude');
%     title('Signal of Interest with Gaussian filter');
% 
%     subplot(3,1,2)
%     plot(b,vg)
%     xlabel('Time [sec]');
%     ylabel('Amplitude');
%     title('Gabor Filtered Signal of Interest');
%     axis([0 9 -1 1])
% 
%     subplot(3,1,3)
%     plot(ks,vgt/max(vgt))
%     xlabel('\omega [Hz]');
%     ylabel('Amplitude');
%     title('Signal of Interest in frequency domain (Gabor)');
%     %axis([-2.6 2.6 0 1])
%     drawnow
    %pause(0.1)
    
end



% figure(4)
% subplot(3,1,1)
% plot(b,v,'b'), hold on
% plot(b,gab,'r','Linewidth',3), hold off
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest with Gaussian filter');
% 
% subplot(3,1,2)
% plot(b,vg)
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Gabor Filtered Signal of Interest');
% axis([0 9 -1 1])
% 
% subplot(3,1,3)
% plot(ks,vgt/max(vgt))
% xlabel('\omega [Hz]');
% ylabel('Amplitude');
% title('Signal of Interest in frequency domain (Gabor)');
%print('3 plot: tau: samples:','-djpeg')


figure(1) = figure('Position', [200, 100, 1250, 1250]);
subplot(1,3,1)
pcolor(tslide,ks/2*pi,(vgt_spec.')), shading interp
title({'Spectrogram of Handle w/ Guassian filter ','(\tau = 100)'});
xlabel('Time [sec]');
ylabel('\omega [Hz]');
set(gca,'Fontsize',[12])
%set(gca,'Ylim',[-20 20])
colormap(hot)
colorbar

subplot(1,3,2)
pcolor(tslide,ks/2*pi,(signal_mexihat_spec.')), shading interp
title({'Spectrogram of Handle w/ Mexican Hat filter ','(\tau = 100)'});
xlabel('Time [sec]');
ylabel('\omega [Hz]');
set(gca,'Fontsize',[12])
%set(gca,'Ylim',[-20 20])
colormap(hot)
colorbar


subplot(1,3,3)
pcolor(tslide,ks/2*pi,(signal_shannon_spec.')), shading interp
title({'Spectrogram of Handle w/ Shannon filter ','(\width = 1 sec)'});
xlabel('Time [sec]');
ylabel('\omega [Hz]');
set(gca,'Fontsize',[12])
%set(gca,'Ylim',[-20 20])
colormap(hot)
colorbar
saveas(gcf,'spectrogram compare: tau = 100, samp step: 0.08925 sec shift.png')


%% Mexican hat filtering for loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
tic

load handel
v = y';

% p8 = audioplayer(v,Fs);
% playblocking(p8);
 

L = (1/Fs)*length(v); % Seconds: 1/(sample rate) * # of samples
n = length(v);
t2 = linspace(0,L,n+1); 
t = t2(1:n);
b = (1:length(v))/Fs;
tau = 100; % Width of Gabor/Gaussian filter

k = (2*pi/(L))*[0:(n/2) -n/2:-1]; ks = fftshift(k); % Frequency domain
signal_mexihat_spec = [];
mexican_hat = (1-(t-4.5).^2).*exp((-tau*(t-4.5).^2)); % Mexican Hat filter
filter_sig = mexican_hat.*v;
filter_sig_fft = abs(fftshift(fft(filter_sig)));


figure(1) = figure('Position',[100, 100, 1250, 1250]);
subplot(3,5,[1 2])
plot(b,v,'b'), hold on
plot(b,mexican_hat,'r','Linewidth',3), hold off
xlabel('Time [sec]');
ylabel('Amplitude');
title('Audio Signal of Interest w/ Mexican hat filter (\tau = 100)');
set(gca,'Fontsize',[12])

subplot(3,5,[6 7])
plot(b,filter_sig)
xlabel('Time [sec]');
ylabel('Amplitude');
title('Mexican Hat Filtered Signal of Interest');
axis([0 9 -1 1])
set(gca,'Fontsize',[12])

subplot(3,5,[11 12])
plot(ks,filter_sig_fft/max(filter_sig_fft))
xlabel('\omega [Hz]');
ylabel('Amplitude');
title('Signal of Interest in Frequency Domain (Mexican hat filtered)');
set(gca,'Fontsize',[12])

tslide = 0:.01:L;
for j = 1:length(tslide)
    
    mexican_hat = (1-(t-tslide(j)).^2).*exp((-tau*(t-tslide(j)).^2)); % Mexican Hat filter
    
    filter_sig = mexican_hat.*v;
    filter_sig_fft = abs(fftshift(fft(filter_sig)));
    signal_mexihat_spec = [signal_mexihat_spec; filter_sig_fft/max(filter_sig_fft)];

    %%% Mexican Hat filtering %%%% Spicy! 
    % Very cool animation of sweeping through signal with filter!
    %figure(5)
    
%     subplot(3,1,1)
%     plot(b,v,'b'), hold on
%     plot(b,mexican_hat,'r','Linewidth',3), hold off
%     xlabel('Time [sec]');
%     ylabel('Amplitude');
%     title('Signal of Interest w/ Mexican hat filter');
% 
%     subplot(3,1,2)
%     plot(b,filter_sig)
%     xlabel('Time [sec]');
%     ylabel('Amplitude');
%     title('Mexican Filtered Signal of Interest');
%     axis([0 9 -1 1])
% 
%     subplot(3,1,3)
%     plot(ks,filter_sig_fft/max(filter_sig_fft))
%     xlabel('\omega [Hz]');
%     ylabel('Amplitude');
%     title('Signal of Interest in Frequency Domain (Mexican hat)');
%     %axis([-2.6 2.6 0 1])
%     drawnow
%     pause(0.01)
end


figure(1)
subplot(3,5,[3 4 5 8 9 10 13 14 15])
pcolor(tslide,ks,(signal_mexihat_spec.')), shading interp
title('Spectrogram of Handle w/ Mexican Hat filter');
xlabel('Time [sec]');
ylabel('\omega [Hz]');
set(gca,'Fontsize',[12])
%set(gca,'Ylim',[-20 20])
colormap(hot)
colorbar
saveas(gcf,'Mex hat: tau = 100, samp step: 0.08925 sec shift.png')



%% Part 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
%  figure(3)
%  subplot(3,1,1)
[y,Fs] = audioread('music1.wav');

L = (length(y)/Fs);
D = round(length(y)/Fs); % Seconds: 1/(sample rate) * # of samples
y1 = downsample(y,5);

piano = y1';
n = length(y1);
x2 = linspace(0,L,n+1); 
x = x2(1:n);

k2 = (2*pi/D)*[0:(n/2)-1 -n/2:-1]; ks2 = fftshift(k2); % Frequency domain for piano

tau = 100; % Width of Gabor/Gaussian filter


tslide = 0:0.1:D;

piano_notes = zeros(length(tslide),1); %Initilizing piano notes vector
piano_gabor2_spec = zeros(length(tslide),length(piano));

for l = 1:length(tslide)

    gabor2 = exp(-tau*(x-tslide(l)).^2); % Gabor/Gaussian filter
    piano_gabor2 = (gabor2.*piano); %
    piano_gabor2_lp = lowpass(piano_gabor2, 370, Fs);
    piano_gabor2_ft = fft(piano_gabor2_lp); % Transformed gabor signal
    piano_gabor2_spec(l,:) = abs(fftshift(piano_gabor2_ft));
    
    [MaxV,Ind] = max(piano_gabor2_ft); %Index where frequency is maximum
    piano_notes(l,:) = abs(k2(Ind))/(2*pi); %maximum index in frequency.
    
    
end


%%% Spectrogram of piano
figure(9)
subplot(2,1,1)
pcolor(tslide, abs(ks2)/(2*pi), piano_gabor2_spec.'), shading interp
title('Spectrogram of piano with lowpass filter (\alpha = 100)');
xlabel('Time [sec]');
ylabel('Frequency [Hz]');
set(gca,'Fontsize',[12])
set(gca,'Ylim',[0 500])
colormap(hot)
colorbar
%saveas(gcf,'Piano spectrogram \alpha = 100.png')



 

[y,Fs] = audioread('music2.wav');

L = (length(y)/Fs);
D = round(length(y)/Fs); % Seconds: 1/(sample rate) * # of samples
y1 = downsample(y,5);

recorder = y1';
n = length(y1);
x2 = linspace(0,L,n+1); 
x = x2(1:n);

k2 = (2*pi/D)*[0:(n/2) -n/2:-1]; ks2 = fftshift(k2); % Frequency domain for piano

tau = 100; % Width of Gabor/Gaussian filter

tslide2 = 0:0.1:D;

recorder_notes = zeros(length(tslide2),1); %Initilizing piano notes vector
recorder_gabor2_spec = zeros(length(tslide2),length(recorder));

for l = 1:length(tslide2)

    gabor2 = exp(-tau*(x-tslide2(l)).^2); % Gabor/Gaussian filter
    recorder_gabor2 = (gabor2.*recorder); %

    recorder_gabor2_ft = fft(recorder_gabor2); % Transformed gabor signal
    [MaxV,Ind] = max(recorder_gabor2_ft); %Index where frequency is maximum
    recorder_notes(l,:) = abs(k2(Ind))/(2*pi); %maximum index in frequency.
    
    recorder_gabor2_spec(l,:) = abs(fftshift(recorder_gabor2_ft));
    
end

% Spectrogram of recorder
figure(9)
subplot(2,1,2)
pcolor(tslide2,abs(ks2)/(2*pi),(recorder_gabor2_spec.')), shading interp
title('Spectrogram of recorder (\alpha = 100)');
xlabel('Time [sec]');
ylabel('frequency [Hz]');
set(gca,'Fontsize',[12])
set(gca,'Ylim',[600 1300])
colormap(hot)
colorbar
saveas(gcf,'Piano & recorder spectrograms.png')
 
 
%% Scores
figure(10) = figure('Position', [100, 100, 900, 900]);
subplot(2,1,1)
plot(tslide,piano_notes,'o','MarkerFaceColor', 'k'); 
yticks([246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23]); 
yticklabels({'B3','C4','C#4','D4','E4','F4'});
ylim ([200 350])
set(gca,'Fontsize',[15])
title("Score for Piano");
grid on
xlabel("Time (s)"); ylabel("Notes");

figure(10)
subplot(2,1,2)
plot(tslide2,recorder_notes,'o', 'MarkerFaceColor', 'k') 
yticks([698.46, 739.99, 783.99, 830.61, 880, 932.3275]); 
yticklabels({'F5','F#5','G5','G#5','A5','A#5'});
ylim ([600 950])
set(gca,'Fontsize',[15])
title("Score for Recorder");
grid on
xlabel("Time (s)"); ylabel("Notes");
saveas(gcf,'Piano & recorder Score.png')
