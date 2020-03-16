function [FileNameList] = Get_audio_CSV_samples()
%%%%%%   AMATH 482   %%%%%%%%%%%%%%%
%%%%%%%   Winter 2020  %%%%%%%%%%%%%
%%%%%%%%   Homework 4    %%%%%%%%%%%  
%%%%%%%%%   Song ID & SVD %%%%%%%%%%
%
% Objective: 
%   1. Process data  to "fit" into the SVD compatible format
%   2. 
%   3. 
%   4. Produce pretty plots and visual data analysis
%   5. 
%
% Procedure:
%   1. load data and inquire about the structure.
%   2. Use structure to get each frame of data. Note: uint8 data type and
%   RGB. [n m RGB frame]
%   3. 
%   4. 
%   5. 
%
% Notes: The following code could be expanded to include more automation
% for...
%   1. Functions for graphing.
% SUMMARY: The following code is a demonstration of .... 

tic
clear all; close all; clc;
cd('/Users/josefhenthorn/Desktop/HW 4/Music folder')
old_folder = pwd;



% %[a,Fs] = webread('url')
% size(a)
% tr_piano=length(a)/Fs; % record time in seconds 
% L = round(length(y)/Fs);
% a = downsample(a,5);




%%%% From audio file

mySelectedFolder = uigetdir(pwd, 'Select a folder'); %gets directory
genre_selector = menu('Enter genre for folder','Indie Rock','Soul','Classic Rock','R&B Hip-hop','Old Folk','Old School Country','Electronic','Modern Pop');
if genre_selector == 1
    genre = "Indie Rock";
elseif genre_selector == 2
    genre = "Soul";
elseif genre_selector == 3
    genre = "Classic Rock";
elseif genre_selector == 4
    genre = "R&B Hip-hop";
elseif genre_selector == 5
    genre = "Old Folk";
elseif genre_selector == 6
    genre = "Old School Country";
elseif genre_selector == 7
    genre = "Electronic";
elseif genre_selector == 8
    genre = "Modern Pop";
end


cd(mySelectedFolder); % Change directory to path of the folder chosen
music_Files_mp3 = dir(fullfile(mySelectedFolder, '*.mp3')); % Load all the music files of this folder.
music_Files_m4a = dir(fullfile(mySelectedFolder, '*.m4a')); % Load all the csv files.


music_Files = vertcat(music_Files_mp3,music_Files_m4a);
filepath = mySelectedFolder;
folder_path = fileparts(mySelectedFolder);


file_names = cell(size(music_Files));
FileNameList = string(zeros(size(music_Files)));

for j = 1:(length(music_Files))
    file_names{j} = {music_Files(j).name};
    
    
    
    
    FileName = string(cell2mat(file_names{j}));
    FileNameList(j) =  string(FileName);
    
    
    file_name = string(string(filepath)+'/'+string(FileName));
    filename = file_name;
    audio_info = audioinfo(file_name);
  
    %j = 1:length(music_Files);
     Artist = audio_info.Artist;
     song_name = audio_info.Title;
     %Artist(j) = getfield(audio_info(j),'Artist'); 
    [y,Fs] = audioread(file_name);
    %size(y)
    L = round(length(y)/Fs); % record time in seconds
    y = downsample(y,8);
    jj = (1:3)*6; % Sample intervals in seconds

    for i = 1:L       
           if i ==  jj(1) &&  i*Fs-2*Fs > 0
               audio_signal1 = y(i*Fs-2*Fs:i*Fs+3*Fs); % Five second sample 1
           end
           if i == jj(2) && i*Fs < length(y)

                audio_signal2 = y(i*Fs-7*Fs:i*Fs-2*Fs); % Five second sample 2
           end   
           if i == jj(3) && i*Fs < length(y)
               audio_signal3 = y(i*Fs-7*Fs:i*Fs-2*Fs); % Five second sample 3
           end
    end
    n = length(audio_signal1);
    LL = round(length(audio_signal1)/Fs);
    x2 = linspace(0,L,n+1); 
    x = x2(1:n);
    k = (2*pi/LL)*[0:(n/2) -n/2:-1]; ks = fftshift(k); % Frequency domain 
    
    samp_step = Fs/2;
    tslide = 0:samp_step:length(audio_signal1);
    t = 1:length(tslide);

    if exist('audio_signal1','var')==1
    Audio_spectrogram1 = zeros(length(tslide),length(ks));
    end
    if exist('audio_signal2','var')==1
    Audio_spectrogram2 = zeros(length(tslide),length(ks));
    end
    if exist('audio_signal3','var')==1
    Audio_spectrogram3 = zeros(length(tslide),length(ks));
    end

    if exist('audio_signal1','var') == 1

        Audio_spectrogram1 = abs(spectrogram(audio_signal1));
        Audio_spectrogram1 = reshape(Audio_spectrogram1,numel(Audio_spectrogram1),1);
    end
    if exist('audio_signal2','var') == 1
            Audio_spectrogram2 = abs(spectrogram(audio_signal2));
            Audio_spectrogram2 = reshape(Audio_spectrogram2,numel(Audio_spectrogram2),1);

    end
    if exist('audio_signal3','var') == 1
        Audio_spectrogram3 = abs(spectrogram(audio_signal3));
        Audio_spectrogram3 = reshape(Audio_spectrogram3,numel(Audio_spectrogram3),1);
        %j
        %Audio_spectrogram_reporter = Audio_spectrogram3(50:75)
    end
    
    folder_path_genre_csv = (string(old_folder)+'/'+genre+' csv data');
    cd(folder_path_genre_csv); % Change directory to path of the folder

     if exist('audio_signal1','var')==1
        writematrix(Audio_spectrogram1,'Sgram_1_'+string(genre)+'_'+string(Artist)+'_'+string(song_name)+'.csv')
     end
     if exist('audio_signal2','var')==1
        writematrix(Audio_spectrogram2,'Sgram_2_'+string(genre)+'_'+string(Artist)+'_'+string(song_name)+'.csv') 
     end
     if exist('audio_signal3','var')==1
        writematrix(Audio_spectrogram3,'Sgram_3_'+string(genre)+'_'+string(Artist)+'_'+string(song_name)+'.csv')
     end
     if Audio_spectrogram1 == Audio_spectrogram2
         j
         flag = "Flag 1: 1&2 R="
     end
end
cd(old_folder)
time = toc
end

