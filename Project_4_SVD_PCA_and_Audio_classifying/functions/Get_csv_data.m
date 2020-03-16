function [genre csv_matrix,csv_content_vect] = Get_csv_data(csv_data_to_load)
% This function loads all csv files within a folder and concatenates them
% into one large column vector.
%

close all; clc
tic
cd('/Users/josefhenthorn/Desktop/HW 4/Music folder')
old_folder = pwd;

csv_data_to_load

myFolder = uigetdir(pwd, 'Select a folder'); %gets directory
cd(myFolder); % Change directory to path of the folder
csvFiles = dir(fullfile(myFolder, '*.csv')); % Load all the csv files.
Number_of_csv = length(csvFiles);

directory = strcat(old_folder,'/');
newStr = cell2mat(split(string(myFolder),string(directory)));
genre = cell2mat(split(newStr," csv data"));

%myFiles = ('*.csv'); %gets all csv files in struct

 %csv_content = readtable('myfile.csv')
for k = 1:length(csvFiles)
    csvFilesnames(k,1) = string(csvFiles(k).name);
end

csv_content_vect = [];

for j = 1:numel(csvFilesnames)

    csv_content = {csvread(csvFilesnames(j,:))};
    csv_content_vect = horzcat(csv_content_vect, csv_content);
end

for jj = 1:numel(csv_content_vect)
    csv_sizes(jj) = size(csv_content_vect{jj},1); 
end

trunc_limit = min(csv_sizes);
csv_matrix = [];

for kk = 1:length(csv_content_vect)
    csv_column = csv_content_vect{kk};
    csv_column = csv_column(1:trunc_limit);
    csv_matrix(:,kk) = csv_column;
end






cd(old_folder)
end


