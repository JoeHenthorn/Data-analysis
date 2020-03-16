%%%%%%%%%      3D Fourier Transform  %%%%%%%%%
%%%%%%%%%         Josef Henthorn     %%%%%%%%%
%%%%%%%%%          Winter 2020       %%%%%%%%% 
%%%%%%%%%   University of Washington %%%%%%%%%
%
% Objective: 
%   1. Multidimensional fft of spatial data
%   2. Frequency filtering
%   3. 3D Fourier analysis
%   4. Pretty plots and gif animation
%
% Procedure:
%   1. load data and make matrices (meshgrid, linspace, reshape, etc.)
%   2. FFT data and avg. then find frequency coordinates of interest.
%   3. Build filter for noisy data with frequency coordinates of interest.
%   4. FFT data again and use filter.
%   5. iFFT data and visualize. 
%
% SUMMARY: The following code is a demonstration of .... 
% 3D spatial Fourier Transform. 
% Clean up noisy signal w/ 3D Gaussian filter.
% Utilize the reshape function to structure a matrix.
% Utilize the ind2sub to find specific values in a matrix.
% Two methods for adding carriage returns to titles.
% Pretty 3D plots (4D including time).
% 3D plotting with a for-loop to simulate time steps (animation).
% Saving figures as pictures
% Capture plots as images and save them as a .gif file.

clear; close all; clc;
tic
load Testdata
filename = 'testAnimated.gif';

L = 15; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x=x2(1:n); y=x; z=x; % three identical linspace vectors
k=(2*pi/(L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); % 3D special frequency-space matrix


tau = 0.2; % Width of Gaussian filter
V = zeros(n,n,n); % Preload a matrix for efficiency

% This for-loop unpacks and structures data (also plots).
% This for-loop performed an FFT to get 3D frequency coordinates of an object.
% This for-loop sums the FFT for each time point.
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   %close all, isosurface(X,Y,Z,abs(Un),0.4)
   %axis([-20 20 -20 20 -20 20]), grid on
   F = (fftn(Un));
   V = (V+F);
  
end

V = abs(V/20); % Average FFT's and use only real values.

[Value, indx] = max(V(:)); % Find the max value (highest amplitude freq.) in the matrix. 
% ind2sub is a coordinates finder.
% Input: matrix size and index value of interest. Output is the coordinates of interest for that value.
[sx, sy, sz] = ind2sub([64,64,64],indx); 

%%% Specific "frequency coordinates" for the object. %%%
%%% These are the k_0 values for the Gaussian filter.
k_0x = Kx(sx, sy, sz); % Plug in specific point coord. into freq.space matrix 
k_0y = Ky(sx, sy, sz); 
k_0z = Kz(sx, sy, sz);


% Gaussian filter: exp(-tau*(k - k_0)^2), tau is width, k0 is the specific
% frequency coordinates and k is the special frequency-space vector 
% (see above)

% Gaussian filter: exp(-tau*( (kx - k_0x)^2 + ... + (kz - k_0z)^2 )
filter = exp( -tau*((Kx - k_0x).^2 + (Ky - k_0y).^2 + (Kz - k_0z).^2));


%%% For pretty plot animation
duration = 20; % Time end
t = 0:1:duration; % Time vector
%%% Preset matrices
Object_position = zeros(1,3);
p = zeros(1,3);


% This for-loop unpacks and structures data (also pretty plots w/ Gif).
% This for-loop performed an FFT and filters the noisy signal w/ Gaussian filter and then performs an INV FFT.
% This for-loop uses clean data to find the object in space and time.
% This for-loop builds up images and makes a .gif file.
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   F2 = (fftn(Un));
   filter_freq = F2.*filter;
   Clean_sig = real(ifftn(filter_freq));
   
   [val, index] = max(Clean_sig(:));
   [q,r,s] = ind2sub([n,n,n],index); % Input: matrix size and index of interest. Output is the coordinates of interest.

   % Spatial positions/ coordinates of object/marble
   X_p = X(q, r, s);
   Y_p = Y(q, r, s);
   Z_p = Z(q, r, s);
   
   Object_position(j,:) = [X_p,Y_p,Z_p]; % Space & time object coordinates matrix. Every row is a new time point of 3D(x,y,z) row vectors.
   
   f = figure(1);
   figure(1)
   p1 = plot3(Object_position(:,1),Object_position(:,2),Object_position(:,3),'b','Linewidth',3); grid on, hold on
   p2 = plot3(Object_position(j,1),Object_position(j,2),Object_position(j,3),'r.','markersize',20);  
   title(sprintf('Tracking an Object in Time And 3D Space \n t = %g', t(j))); % Timing counting title
   xlabel('X axis');
   ylabel('Y axis');
   zlabel('Z axis');
   axis([-11 16 -6 6 -9 11]);
   p(j,:) = [Object_position(j,1) Object_position(j,2) Object_position(j,3)];
   location = sprintf('Marble Location: \n x: %.2g \n y: %.2g \n z: %.2g',p(j,:));
   legend('Marble trajectory',location)
   set(gca,'Fontsize',[12])
   drawnow
   pause(0.0001)
   
   % Capture the plot as an image 
   frame = getframe(f); 
   im = frame2im(frame); 
   [imind,cm] = rgb2ind(im,256); 
   % Write to the GIF File 
   if j == 1 
       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
   else 
       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
   end 
   hold off
  
end

% Final object position 3D pretty plot 
figure(1)
plot3(Object_position(:,1),Object_position(:,2),Object_position(:,3),'b','Linewidth',3), grid on, hold on
plot3(Object_position(20,1),Object_position(20,2),Object_position(20,3),'r.','markersize',20), % Final location of object
title({'Tracking an Object in Time And 3D Space','t = 20'})
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
p = [Object_position(20,1) Object_position(20,2) Object_position(20,3)];
location = sprintf('Marble Location: \n x: %.2g \n y: %.2g \n z: %.2g',p);
legend('Marble trajectory',location)
axis([-11 16 -6 6 -9 11]);
set(gca,'Fontsize',12)
print('Marble trajectory & location','-djpeg')

Last_point = [Object_position(20,1) Object_position(20,2) Object_position(20,3)];


% % % % % % %% The following plots are for visualizing the signal, Fourier analysis, and filter in 3D
% % % figure(2) % Noisy ultrasound data in space
% % % isosurface(X,Y,Z,abs(Un),0.4)
% % % %%
% % % figure(3) % Noisy ultrasound data in frequency
% % % 
% % % isosurface(Kx,Ky,Kz,(abs(V)/real(max(V))),0.6), hold on
% % % 
% % % %%
% % % figure(4) % "Clean" ultrasound data in space
% % % isosurface(X,Y,Z,Clean_sig,0.6)
% % % %%
% % % figure(5) % "Clean" ultrasound data in frequency
% % % isosurface(Kx,Ky,Kz,abs(filter_freq)./max(filter_freq(:)),1)
% % % %%
% % % figure(8) % 3D Guassian Filter
% % % isosurface(Kx,Ky,Kz,abs(filter),0.6)
% % % 
% % % %%



%% The following code demonstrates frequency filtering for each frequency "direction."

 
figure(6) % Frequency space averaged signals 
Vx = V(1,1,:);       % 64     1
Vx = squeeze(Vx)/max(Vx);
Vy = V(1,:,1)';       % 1     1    64
Vy = squeeze(Vy)/max(Vy);    % 64     1
Vz = V(:,1,1);      % 64     1
Vz = squeeze(Vz)/max(Vz);


% Gaussian filter EA direction: exp(-tau*( (kx - k_0x)^2 + ... + (kz - k_0z)^2 )
filter_x = exp( -tau*((ks - k_0x).^2 ));

filter_y = exp( -tau*((ks - k_0y).^2 ));

filter_z = exp( -tau*((ks - k_0z).^2));


b = linspace(-L,L,n);
t = linspace(0,20,20);


subplot(3,1,1)
plot(b,fftshift(Vx),'r','Linewidth',2), hold on, grid on
title('Signal Frequencies in \omega_x with Guassian filter') 
subplot(3,1,1), hold on
plot(b,filter_y,'k','Linewidth',2)
legend('Noisy signal x','filter','Location','Best')
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)

subplot(3,1,2), hold on
plot(b,fftshift(Vy),'g','Linewidth',2), grid on
title('Signal Frequencies in \omega_y with Guassian filter') 
subplot(3,1,2)
plot(b,filter_x,'k','Linewidth',2), hold on, grid on
legend('Noisy signal y','filter')
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)
   
subplot(3,1,3)
plot(b,fftshift(Vz),'b','Linewidth',2), hold on, grid on
title('Signal Frequencies in \omega_z with Guassian filter') 
subplot(3,1,3)
plot(b,filter_z,'k','Linewidth',2), hold on, grid on, grid on
legend('Noisy signal z','filter')
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)
print('Noisy frequency + filters','-djpeg')


figure(7) % Frequency filtered  signals 
FF = abs(filter_freq);
FFx = FF(:,1,1);       % 64     1
FFy = FF(1,1,:);       % 1     1    64
FFy = squeeze(FFy);    % 64     1
FFz = FF(1,:,1)';      % 64     1
FFz = squeeze(FFz);



subplot(3,1,1)
plot(b,(FFx/max(FFx)),'r','Linewidth',2), hold on, grid on
title('Filtered signal in \omega_x Domain') 
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)
   
subplot(3,1,2)
plot(b,FFz/max(FFz),'g','Linewidth',2), hold on, grid on
title('Filtered signal in \omega_y Domain') 
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)
   
subplot(3,1,3)
plot(b,FFy/max(FFy),'b','Linewidth',2), hold on, grid on
title('Filtered signal in \omega_z Domain') 
xlabel('frequency');
ylabel('amplitude');
set(gca,'Fontsize',12)
print('Filtered frequency signal','-djpeg')

%% Filtered Spatial data plots 

figure(8)
subplot(3,1,1)
spatial1 = fftshift(abs(ifft(FFx)));
plot(b,spatial1,'r','Linewidth',2), hold on, grid on
title('Filtered signal in spatial domain') 
xlabel('x');
ylabel('magnitude');
set(gca,'Fontsize',12)
   
subplot(3,1,2)
spatial2 = fftshift(abs(ifft(FFy)));
plot(b,spatial2,'g','Linewidth',2), hold on, grid on
%title('Filtered signal in spatial domain') 
xlabel('y');
ylabel('magnitude');
set(gca,'Fontsize',12)
   
subplot(3,1,3)
Spatial3 = fftshift(abs(ifft(FFz)))/L;
plot(b,Spatial3,'b','Linewidth',2), hold on, grid on
%title('Filtered signal in spatial domain') 
xlabel('z');
ylabel('magnitude');
set(gca,'Fontsize',12)
print('Filtered spatial signal','-djpeg')

%%

figure(9)
filt = exp( -0.2*(k).^2);
plot(b,fftshift(filt),'k','Linewidth',2), grid on 
title('Gaussian Curve') 
xlabel('\omega');
ylabel('amplitude');
axis([-15 15 0 1.1])
set(gca,'Fontsize',12)
print('Gaussain Filter','-djpeg')

run_time = toc + "  sec."