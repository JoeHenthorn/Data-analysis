%%%%%%%%%        Neural Networks     %%%%%%%%%
%%%%%%%%%         Josef Henthorn     %%%%%%%%%
%%%%%%%%%          Winter 2020       %%%%%%%%% 
%%%%%%%%%   University of Washington %%%%%%%%%

% Objective: 
%   1. load in data set 
%   2. Convert data and make training and validation groups
%   3. Build a fully connected network and observe results
%   4. Build a convolutional network and observe results
%   5. Produce pretty plots and visual data analysis
%
% SUMMARY: The following code is a demonstration of .... 
%  Fully connected and convultional neural networks.
% Be sure to check the Deep Network Designer in MatLab for a GUI network
% builder with many visual examples. 

clear all; close all; clc;
tic
cd('/Users/mycomputername/Desktop/Neural_Networks')
old_folder = pwd;


load('fashion_mnist.mat')


%tic
[m n p] = size(X_train);
% [X_test] = UTILITY_wavelet(X_train);
% time_wavelet = toc

% 
% figure(1) % Show the first items in the data set
% for k = 1:9
%     subplot(3,3,k)
%     imshow(reshape(X_train(k,:,:),[n p]));
% end



 
%%

X_train = im2double(X_train);
X_train = reshape(X_train,60000,28, 28,1);
X_train = permute(X_train,[2 3 4 1]);


X_test = im2double(X_test);
X_test = reshape(X_test,10000,28,28,1);

X_test = permute(X_test,[2 3 4 1]);


X_validate = X_train(:,:,:,1:5000);
X_train = X_train(:,:,:,5001:end);




Y_validate = categorical(y_train(1:5000))';
y_train = categorical(y_train(5001:end))';
y_test = categorical(y_test)';

Epochs = 5;

%% %%%%%% Building a fully-connected Neural Network 
Epochs = 10; % Number of cycles




first_layer = 555; % Number neurons in 1st layer
second_layer = 200; % Number neurons in 2nd layer
third_layer = 300; % Number neurons in 3rd layer
fourth_layer = 100; % Number neurons in 4th layer
fifth_layer = 100; % Number neurons in 5th layer
sixth_layer = 10; % Number neurons in 6th layer


% layers and the activation functions.
layers = [imageInputLayer([n p 1])
        fullyConnectedLayer(first_layer)
        reluLayer
%         fullyConnectedLayer(second_layer)
%         reluLayer
%         fullyConnectedLayer(third_layer)
%         reluLayer
%         fullyConnectedLayer(fourth_layer)
%         reluLayer
%         fullyConnectedLayer(fifth_layer)
%         reluLayer
%         fullyConnectedLayer(sixth_layer)
%         reluLayer
        fullyConnectedLayer(10)
        softmaxLayer
        classificationLayer];

   
% Specify hyperparameters of the neural network.
options = trainingOptions('adam', ... % Creg Gin recommended optimizer
    'MaxEpochs',Epochs,... % "Rounds" of training 
    'InitialLearnRate',1e-3, ...% step size in gradient decent ??
    'L2Regularization',1e-4, ...% Two norm 
    'ValidationData',{X_validate,Y_validate}, ...%  validation
    'Verbose',false, ... % Don't print to screen
    'Plots','training-progress'); % Plot progress


% Train the system.
network = trainNetwork(X_train,y_train,layers,options);


%currentfig = findall(groot,'Type','Figure');
%savefig(gcf,'Training Progress: Epoch = 5: Network 100_100_.png')

%% Confusion for training
figure(2)
y_pred = classify(network,X_validate);

plotconfusion(Y_validate,y_pred)
title('F.C. NN Confusion Matrix training')
saveas(gcf,'Training F.C Confusion Matrix : Epoch'+string(Epochs)+' Network.png')

y_pred2 = classify(network,X_train);
figure(3)
plotconfusion(y_train,y_pred2)
title('F.C. NN Confusion Matrix validation')
saveas(gcf,'Validation F.C Confusion Matrix : Epoch'+string(Epochs)+' Network.png')

% Test classifier
y_pred3 = classify(network,X_test);
figure(4)
plotconfusion(y_test,y_pred3)
title('F.C. NN Confusion Matrix Test data')
saveas(gcf,'Test F.C Confusion Matrix : Epoch'+string(Epochs)+' Network.png')




%% %%%%% Convolutional Neural network. CNN   %%%%%%%%%%%%%%%%%%%%%%%%%
Epochs = 5;

layers = [
    imageInputLayer([28 28 1],"Name","imageinput")
    convolution2dLayer([5 5],100,"Name","conv_1","Padding","same")
    tanhLayer("Name","tanh_1")
    averagePooling2dLayer([2 2],"Name","avgpool2d_1","Padding","same","Stride",[2 2])
    convolution2dLayer([5 5],200,"Name","conv_2")
    tanhLayer("Name","tanh_2")
    averagePooling2dLayer([2 2],"Name","avgpool2d_2","Padding","same","Stride",[2 2])
    convolution2dLayer([5 5],300,"Name","conv_3")
    reluLayer
    fullyConnectedLayer(100,"Name","fc_1")
    reluLayer
    fullyConnectedLayer(10,"Name","fc_2")
    softmaxLayer("Name","softmax")
    classificationLayer("Name","clxassoutput")];


% Specify hyperparameters of the neural network.
options = trainingOptions('adam', ... % Creg Gin recommended optimizer
    'MaxEpochs',Epochs,... % "Rounds" of training 
    'InitialLearnRate',1e-3, ...% step size in gradient decent ??
    'L2Regularization',1e-4, ...% Two norm 
    'ValidationData',{X_validate,Y_validate}, ...%  validation
    'Verbose',false);%, ... % Don't print to screen
    %'Plots','training-progress'); % Plot progress


% Train the system.
network_CNN = trainNetwork(X_train,y_train,layers,options);

% figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
% plot(layerGraph(layers))
%%

%% Confusion LeNet 5
figure(6)

y_pred = classify(network_CNN,X_validate);
plotconfusion(Y_validate,y_pred)
title('Custom CNN Confusion Matrix training')
%saveas(gcf,'Training Custom CNN Confusion Matrix : Epoch'+string(Epochs)+' Network'+string(second_layer)+'.png')
figure(7)
y_pred = classify(network_CNN,X_train);
plotconfusion(y_train,y_pred)
title('Custom Confusion Matrix validation')
%saveas(gcf,'Validation Custom CNN Confusion Matrix validation: Epoch'+string(Epochs)+' Network'+string(second_layer)+'.png')
% Test classifier
figure(8)
y_pred = classify(network_CNN,X_test);
plotconfusion(y_test,y_pred)
title('Custom Confusion Matrix Test data')
%saveas(gcf,'Test Custom CNN Confusion Matrix: Epoch'+string(Epochs)+' Network'+string(second_layer)+'.png')




%%
% 
% 
% 
% % Artificial Neurons aka...
% layer_nodes = [5];
% 
% % Initiate the Feedforward network
% network = feedforwardnet(layer_nodes);
% 
% % Feed in data for training
% network = train(network,X_train,y_train);
% 
% 
% 
% %%
% layers1_nodes5_NN = network;
% save layers1_nodes5_NN
% 
% %ypredict = network(x2);
% 
% %plot(x2,ypredict,':r','Linewidth',3)
% %legend('data','Function','Neural Network')
% 
% 
% %load Trainned_NN
% %newoutput = Trainned_NN(newinput);
% 
% 
% 
