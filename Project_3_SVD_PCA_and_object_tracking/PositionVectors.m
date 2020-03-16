function [Position_x,Position_y] = PositionVectors(input)
%No input needed. A GUI will prompt for an input file.
% Output: A vectors of positions for tracked object.  
%
% Objective: 
%   1. Track an object in time and space via video data.
%   2. Return object position data

%
% Procedure:
%   1. Load data and make matrices (meshgrid, linspace, reshape, etc.)
%   2. Binerize the imaje into a logical array
%   3. Average data and isolate object from the rest of frame.

close all;

%prompt = str2num(cell2mat(inputdlg({'Enter experiment number [1-4]:'})))
file_name = 'cam'+string(input)+'.mat';

% for c = 1:3
% if prompt == 1
%     file_name = 'cam'+string(c)+'_1.mat'
%     load(file_name)
% elseif prompt == 2
%     file_name = 'cam'+string(c)+'_2.mat';
%     load(file_name)
% elseif prompt == 3
%     file_name = 'cam'+string(c)+'_3.mat';
%     load(file_name)
% elseif prompt == 4
%     file_name = 'cam'+string(c)+'_4.mat';
%     load(file_name)  
% end
% end

%prompt = {'Enter vidFramesX_Y:'};
%[file_name] = uigetfile('*.mat');
Position_x = []; 
Position_y = []; 


    datastruct = (load(file_name));
    [vidFrames] = v2struct(datastruct);

    n = length(vidFrames(:,1,1,1));
    m = length(vidFrames(1,:,1,1));

    Px = [];
    Py = [];

    numFrames = size(vidFrames,4);
    for j = 1:numFrames
        indx_x = [];
        indx_y = [];

        img1 = rgb2gray((vidFrames(:,:,:,j)));

        %%%%%%%% BELOW IS USED FOR TRUNCATION 
        %%%%%%%% DON'T USE WITH SIDEWAYS VIDS
        img1(:,1:250) = 0; img1(:,450:end) = 0;
        
        img = imbinarize(img1,0.99);
        if img(:,:,:) == 0
            img = imbinarize(img1,0.95);
             if img(:,:,:) == 0
                img = imbinarize(img1,0.90);
                if img(:,:,:) == 0
                    img = imbinarize(img1,0.85);
    %                 if img(:,:,:) == 0
    %                      img = img1;
    %                 end
                end
            end
        end
       % imshow(img)

        for jj = 1:n
            for kk = 1:m
                if img(jj,kk) == true
                    indx_y = [indx_y jj];
                    indx_x = [indx_x kk];
                end
            end
        end
        Px = [Px mean(indx_x)];
        Py = [Py mean(indx_y)];

    end

    P_y = (n - Py);

    %plot(Px,P_y)
    %axis([0 m 0 n])

%    U_input = (inputdlg("Do you want to see the object and tracking? (Y/N)",'s'));
% 
%     if U_input == "y" || U_input == "Y" || U_input == "yes" || U_input == "Yes"|| U_input == " "
% 
%         for j = 1:length(Px) 
% 
%             figure(2)
%             X = rgb2gray(vidFrames(:,:,:,j));
% 
%             imshow(X), hold on
% 
%             plot(Px(j),Py(j),'ro','Markersize',40), hold on
%             title({'Tracking an Object in Time And 3D Space'})
%             xlabel('X axis');
%             ylabel('Y axis');
%             axis([0 m 0 n])
%             set(gca,'Linewidth',[20])
%             drawnow
%         end
%     end

%      Px = Px';
%      Px = Px(~isnan(Px))';
%      
%      Py = Py';
%      Py = Py(~isnan(Py))';
    

    % The code below gets rid of NaN values and replaces them with the
    % value of zero. 
    for j = 1:length(Px) 
        if isnan(Px(j))
            Px(j) = 0;
        end
        if isnan(Py(j))
            Py(j) = 0;
        end
    end
    
    % The code below gets rid of 0 values and replaces them with an
    % interpolated value simular to the nearst value. Beyond this point the motion data is artificialy
    % augmented, and is therefore partially idealized.
    v = (Px == 0);
    Px(v) = interp1(find(~v), Px(~v), find(v), 'nearest','extrap');

    z = (Py == 0);
    Py(z) = interp1(find(~z), Py(~z), find(z), 'nearest','extrap');
    
    
    Position_x = Px-mean(Px,2);
    Position_y= Py-mean(Py,2);
end
    
    

