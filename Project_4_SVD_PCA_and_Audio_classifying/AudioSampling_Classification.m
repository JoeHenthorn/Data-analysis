%%%%%%%%%         Audio sampling      %%%%%%%%%
%%%%%%%%%        &  classification    %%%%%%%%%
%%%%%%%%%         Josef Henthorn      %%%%%%%%%
%%%%%%%%%          Winter 2020        %%%%%%%%% 
%%%%%%%%%   University of Washington  %%%%%%%%%
%
% Objective: 
%   1. Build function to sample audio data
%   2. Make spectrograms with the sample and make csv files of spectrogram data
%   3. Make function to load csv files and match csv data size
%   4. Use the SVD to build a classifier.
%   4. Pretty plots
%
% Procedure:
%   1. load data and make matrices (meshgrid, linspace, reshape, etc.)
%   2. See README 
%
% SUMMARY: The following code is a demonstration of .... 
% Audio Sampling and classification with SVD 


clear all, close all; clc
tic
cd('/Users/mycomputername/Desktop/AudioSampling_Classification')
old_folder = pwd;

%FileNameList = Get_audio_CSV_samples();

% Genre classifier Fill in your own data

load_data = "Classic_Rock csv data";
[genre1,Classic_Rock] = Get_csv_data(load_data);
load_data = "Soul csv data";
[genre2,Soul] = Get_csv_data(load_data);
load_data = "Country csv data";
[genre3,Country] = Get_csv_data(load_data);
load_data = "Electronic csv data";
[genre4,Electronic] = Get_csv_data(load_data);
load_data = "Modern_pop csv data";
[genre5,Modern_pop] = Get_csv_data(load_data);
load_data = "R_&_B csv data";
[genre6,R_and_B] = Get_csv_data(load_data);
load_data = "Indie csv data";
[genre7,Indie] = Get_csv_data(load_data);





%%

%  Finds the three genres with the most data for SVD & makes them the same length

Genre_array = {Classic_Rock, Soul, Country, Electronic, Modern_pop, R_and_B, Indie};

% for j = 1:length(Genre_array)
%     Genre_array(1)
%     
% end





%%%% finding three genres with the most data
len1 = size(Classic_Rock);
len1 = len1(2);
len2 = size(Soul);
len2 = len2(2);
len3 = size(Country);
len3 = len3(2);
len4 = size(Electronic);
len4 = len4(2);
len5 = size(Modern_pop);
len5 = len5(2);
len6 = size(R_and_B);
len6 = len6(2);
len7 = size(Indie);
len7 = len7(2);

length_of_raw_csv = [len1 len2 len3 len4 len5 len6 len7];
[largest3, indexoflargest3] = maxk(length_of_raw_csv, 3);


limit_samp = min(largest3); % Controls the number of files
%limit_samp = 81; % Controls the number of files


for i = 1:length(indexoflargest3)
    
    if indexoflargest3(i) ==1
        if i == 1
            first_genre = "Classic_Rock";
            genre_1 = Classic_Rock;
        elseif i == 2
            second_genre = "Classic_Rock";
            genre_2 = Classic_Rock;
        elseif i == 3
            third_genre = "Classic_Rock";
            genre_2 = Classic_Rock;
        end
    elseif indexoflargest3(i) ==2
        if i == 1
            first_genre = "Soul";
            genre_1 = Soul;
        elseif i == 2
            second_genre = "Soul"; 
            genre_2 = Soul;
        elseif i == 3
            third_genre = "Soul";
            genre_3 = Soul;
        end
    elseif indexoflargest3(i) ==3
        if i == 1
            first_genre = "Country";
            genre_1 = Country;
        elseif i == 2
            second_genre = "Country"; 
            genre_2 = Country;
        elseif i == 3
            third_genre = "Country";
            genre_3 = Country;
        end
    elseif indexoflargest3(i) ==4
        if i == 1
            first_genre = "Electronic";
            genre_1 = Electronic;
        elseif i == 2
            second_genre = "Electronic"; 
            genre_2 = Electronic;
        elseif i == 3
            third_genre = "Electronic";
            genre_3 = Electronic;
        end
    elseif indexoflargest3(i) ==5
        if i == 1
            first_genre = "Modern_pop";
            genre_1 = Modern_pop;
        elseif i == 2
            second_genre = "Modern_pop"; 
            genre_2 = Modern_pop;
        elseif i == 3
            third_genre = "Modern_pop";
            genre_3 = Modern_pop;
        end
    elseif indexoflargest3(i) ==6
        if i == 1
           first_genre = "R_and_B";
           genre_1 = R_and_B;
        elseif i == 2
           second_genre = "R_and_B";
           genre_2 = R_and_B;
        elseif i == 3
           third_genre = "R_and_B";
           genre_3 = R_and_B;
        end
    elseif indexoflargest3(i) ==7
        if i == 1
           first_genre = "Indie";
           genre_1 = Indie;
        elseif i == 2
           second_genre = "Indie"; 
           genre_2 = Indie;
        elseif i == 3
           third_genre = "Indie";
           genre_3 = Indie;
        end
    end
end

Genres_to_Use = [first_genre second_genre third_genre];
genre_1 = genre_1(:,1:limit_samp);
genre_2 = genre_2(:,1:limit_samp);
genre_3 = genre_3(:,1:limit_samp);

%%%% truncating the csv data to match
len_1 = length(genre_1);
len_2 = length(genre_2);
len_3 = length(genre_3);
length_genres = [len_1 len_2 len_3];

limit_spectrogram = min(length_genres);



genre_1 = genre_1(1:limit_spectrogram,:)';
genre_2 = genre_2(1:limit_spectrogram,:)';
genre_3 = genre_3(1:limit_spectrogram,:)';


writematrix(genre_1,'Genre 1.csv')
writematrix(genre_2,'Genre 2.csv')
writematrix(genre_3,'Genre 3.csv')
%%
clear all, close all; clc
cd('/Users/josefhenthorn/Desktop/HW 4/Music folder')
old_folder = pwd;

genre_1 = load('Genre 1.csv');
genre_2 = load('Genre 2.csv');
genre_3 = load('Genre 3.csv');


%%
clear figure

figure(1) = figure('Position', [200, 100, 900, 500]);
subplot(1,3,1)
pcolor(genre_1),shading interp
colorbar
title({'Spectrogram of Genre 1',first_genre});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])

subplot(1,3,2)
pcolor(genre_2), shading interp
colorbar
title({'Spectrogram of Genre 2',second_genre});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])

subplot(1,3,3)
pcolor(genre_3), shading interp
colorbar
title({'Spectrogram of Genre 3',third_genre});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])
saveas(gcf,'Spectrogram of Genre test.png')
clear figure

%%   This is where the magic happens...
for ii = 1:40
num_of_files = size(genre_1);
numfiles = num_of_files(1);
if mod(numfiles,2) ~= 0
    numfiles = numfiles;
end

numtest = 10;


train_set = [1:numfiles];
P = randperm(numfiles,numtest);
P = sort(P);

for i = 1:numtest
    train_set = train_set(train_set~=P(i));
end


train_gen1 = (genre_1(:,train_set));
%size(train_gen1)
train_gen2 = (genre_2(:,train_set));
train_gen3 = ((genre_3(:,train_set)));

test_gen1 = genre_1(:,[P]);
test_gen2 = genre_2(:,[P]);
test_gen3 = genre_3(:,[P]);

feature = numfiles;

[result,w,U,S,V,threshold1,threshold2,sortgen1,sortgen2,sortgen3] = Music_genre_trianer3(train_gen1,train_gen2,train_gen3,feature);

test_set = ([test_gen1 test_gen2 test_gen3]); 

Test_Matrix = U*test_set; % SVD  projection
pval = w'*Test_Matrix; % LDA  projection

testSet = size(test_set);
testSet = testSet(2);

%pvals = sort(pval);
pvals = (pval);

label = zeros(testSet,1);
for i = 1:testSet
    if pvals(1,i) < threshold1
        label(i) = 1;
    elseif pvals(1,i) > threshold1 && pvals(1,i) < threshold2
        label(i) = 2;
    elseif pvals(1,i) > threshold2
        label(i) = 3;
    end
end


hidden = zeros(testSet,1);
hidden(1:(testSet/3)) = 1;
hidden((testSet/3)+1:2*(testSet/3)) = 2;
hidden(2*(testSet/3)+1:testSet) = 3;
total = testSet;


% errNum = sum(abs(label - hidden));
% 
% disp('Rate of success');
% sucRate = 1-errNum/total


correct = zeros(testSet,1);
for i = 1:testSet
    if label(i) == hidden(i)
        correct(i) = 1;
        %indic = [label(i) hidden(i)]
    else
        correct(i) = 0;
    end
end
% close all,
 accuracy(ii) = sum(correct)/total;
 error = 1 - accuracy;
% plot(label,'g'),hold on
% plot(hidden)
end 

Average_accuracy = sum(accuracy)/length(accuracy)


close all,
figure(5)
p1 = plot(sortgen1,0,'ro');hold on
p2 =plot(sortgen2,1,'bo');
p3 =plot(sortgen3,2,'go');

p4 = plot([threshold1 threshold1],[0 3],'k','Linewidth',3);
p5 = plot([threshold2 threshold2],[0 3],'k','Linewidth',3);
title('LDA: Three genre groups');
grid on
xlabel('Projected axis');
ylabel('speration of groups');
set(gca,'Fontsize',[12])
saveas(gcf,'Test 3 genres: Sample data grouping and thresholds.png')



%%
clear figure
figure(6) = figure('Position', [200, 100, 700, 500]);


subplot(1,3,1)
histogram(sortgen1), hold on
title('Projections of test genre 1');
xlabel('Smallest-to-largest projection values');
plot([threshold1 threshold1],[0 30],'r','Linewidth',3);


subplot(1,3,2)
histogram(sortgen2), hold on
title('Projections of test genre 2');
xlabel('Smallest-to-largest projection values');
plot([threshold1 threshold1],[0 30],'r','Linewidth',3);

%plot([threshold2 threshold2],[0 3],'k');

subplot(1,3,3)
histogram(sortgen3), hold on
title('Projections of test genre 3');
xlabel('Smallest-to-largest projection values');
plot([threshold2 threshold2],[0 30],'r','Linewidth',3);
saveas(gcf,'Test 3 genres: Histogram of projected values.png')

%%











