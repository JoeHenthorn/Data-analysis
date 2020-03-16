%%%%%%%%% Oject tracking & Analysis  %%%%%%%%%
%%%%%%%%%           SVD & PCA        %%%%%%%%%
%%%%%%%%%         Josef Henthorn     %%%%%%%%%
%%%%%%%%%          Winter 2020       %%%%%%%%% 
%%%%%%%%%   University of Washington %%%%%%%%% 
%
% Objective: 
%   1. Use a position finding function to track an object
%   2. Process data  to "fit" into the SVD compatible format
%   3. Get principle components and recreate the object motion
%   4. Produce pretty plots and visual data analysis
%
% Procedure:
%   1. load data and inquire about the structure.
%   2. Use structure to get each frame of data. Note: uint8 data type and
%   RGB. [n m RGB frame]
%
% Notes: The following code could be expanded to include more automation
% for...
%   1. Functions for graphing.
% SUMMARY: The following code is a demonstration of .... 

clear all; close all; clc;

%%%%%%%%%%%%%   Test 1  %%%%%%%%%%%%%
close all;

vid = "1_1"
[Posx_1,Posy_1] = PositionVectors("1_1");
% Position vectors for x
% Position vectors for y
vid = "2_1"
[Posx_2,Posy_2] = PositionVectors("2_1");
% Position vectors for x
% Position vectors for y
vid = "3_1"
[Posx_3,Posy_3] = PositionVectors("3_1");
% Position vectors for x
% Position vectors for y
% Position vectors for x & y camera 3
% Reversed for sideways case

Position_Data = {Posx_1;Posy_1;Posx_2;Posy_2;Posy_3;Posx_3};
[singularValues, energy,energy_r_cum,plot11,plot12,plot13] = SVD_Data(Position_Data,1);

singularValues1 = singularValues;
energy1 = energy;
energy1_cum = energy_r_cum';

%%
%%%%%%%%%%%%%   Test 2  %%%%%%%%%%%%%
close all;

vid = "1_2"
[Posx_1,Posy_1] = PositionVectors('1_2');
% Position vectors for x
% Position vectors for y

vid = "2_2"
[Posx_2,Posy_2] = PositionVectors("2_2");
% Position vectors for x
% Position vectors for y

vid = "3_2"
[Posx_3,Posy_3] = PositionVectors("3_2");
% Position vectors for x
% Position vectors for y

% Position vectors for x & y camera 3
% Reversed for sideways case

Position_Data = {Posx_1;Posy_1;Posx_2;Posy_2;Posy_3;Posx_3};

[singularValues, energy,energy_r_cum,plot21,plot22,plot23] = SVD_Data(Position_Data,2);

singularValues2 = singularValues;
energy2 = energy;
energy2_cum = energy_r_cum';

%%
%%%%%%%%%%%%%   Test 3   %%%%%%%%%%%%%
close all;

vid = "1_3"
[Posx_1,Posy_1] = PositionVectors("1_3");
% Position vectors for x
% Position vectors for y

vid = "2_3"
[Posx_2,Posy_2] = PositionVectors("2_3");
% Position vectors for x
% Position vectors for y

vid = "3_3"
[Posx_3,Posy_3] = PositionVectors("3_3");
% Position vectors for x
% Position vectors for y

% Position vectors for x & y camera 3
% Reversed for sideways case

Position_Data = {Posx_1;Posy_1;Posx_2;Posy_2;Posy_3;Posx_3};
[singularValues, energy,energy_r_cum,plot31,plot32,plot33] = SVD_Data(Position_Data,3);

singularValues3 = singularValues;
energy3 = energy;
energy3_cum = energy_r_cum';

%%
%%%%%%%%%%%%%  Test 4  %%%%%%%%%%%%%
close all;

vid = "1_4"
[Posx_1,Posy_1] = PositionVectors("1_4");
% Position vectors for x
% Position vectors for y

vid = "2_4"
[Posx_2,Posy_2] = PositionVectors("2_4");
% Position vectors for x
% Position vectors for y

vid = "3_4"
[Posx_3,Posy_3] = PositionVectors("3_4");
% Position vectors for x
% Position vectors for y

% Position vectors for x & y camera 3
% Reversed for sideways case

Position_Data = {Posx_1;Posy_1;Posx_2;Posy_2;Posy_3;Posx_3};

[singularValues, energy,energy_r_cum, plot41,plot42,plot43] = SVD_Data(Position_Data,4);

singularValues4 = singularValues;
energy4 = energy;
energy4_cum = energy_r_cum';
%%


figure(11)

subplot(1,3,1)

subplot(1,3,1)