
clear all, close all, clc
cd('/Users/josefhenthorn/Desktop/HW 4/Music folder')
old_folder = pwd;



load_data = "Lord Huron";

[Artist1,Lord_Huron] = Get_csv_data(load_data);
cd(old_folder)
load_data = "ODESZA";
[Artist2,ODESZA] = Get_csv_data(load_data);
cd(old_folder)
load_data = "Percy Sledge";
[Artist3,Percy_Sledge] = Get_csv_data(load_data);
cd(old_folder)



%%
len1 = size(Lord_Huron);
len1 = len1(2);
len2 = size(ODESZA);
len2 = len2(2);
len3 = size(Percy_Sledge);
len3 = len3(2);


length_of_raw_csv = [len1 len2 len3];
[largest3, indexoflargest3] = maxk(length_of_raw_csv, 3);


limit_samp = min(largest3); % Controls the number of files
%limit_samp = 81; % Controls the number of files


Artist1 = Artist1(:,1:limit_samp);
Artist2 = Artist2(:,1:limit_samp);
Artist3 = Artist3(:,1:limit_samp);

%%%% truncating the csv data to match
len_1 = length(genre_1);
len_2 = length(genre_2);
len_3 = length(genre_3);
length_genres = [len_1 len_2 len_3];

limit_spectrogram = min(length_genres);

Artist1 = Artist1(1:limit_spectrogram,:)';
Artist2 = Artist2(1:limit_spectrogram,:)';
Artist3 = Artist3(1:limit_spectrogram,:)';

% Save the csv files
writematrix(Artist1,'Test_1_Artist1.csv')
writematrix(Artist2,'Test_1_Artist2.csv')
writematrix(Artist3,'Test_1_Artist3.csv')


%%
clear all, close all; clc
cd('/Users/josefhenthorn/Desktop/HW 4/Music folder')
old_folder = pwd;

Artist1 = load('Test_1_Artist1.csv');
Artist2 = load('Test_1_Artist2.csv');
Artist3 = load('Test_1_Artist3.csv');

%%
figure(1) = figure('Position', [200, 100, 900, 500]); % Remember Clear fig
subplot(1,3,1)
pcolor(Artist1),shading interp
colorbar
title({'Spectrogram of Artist 1','Lord Huron'});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])

subplot(1,3,2)
pcolor(Artist2), shading interp
colorbar
title({'Spectrogram of Artist 2','ODESZA'});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])

subplot(1,3,3)
pcolor(Artist3), shading interp
colorbar
title({'Spectrogram of Artist 3','Percy Sledge'});
xlabel('Time');
ylabel('Frequency');
set(gca,'Fontsize',[12])
saveas(gcf,'Spectrogram of artists test 1.png')
clear figure




%%   This is where the magic happens...
for ii = 1:35
    num_of_files = size(Artist1);
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

    train_gen1 = (Artist1(:,train_set));
    %size(train_gen1)
    train_gen2 = (Artist2(:,train_set));
    train_gen3 = ((Artist3(:,train_set)));

    test_gen1 = Artist1(:,[P]);
    test_gen2 = Artist2(:,[P]);
    test_gen3 = Artist3(:,[P]);

    feature = numfiles;

    [result,w,U,S,V,threshold1,threshold2,sortgen1,sortgen2,sortgen3] = Music_genre_trainer(train_gen1,train_gen2,train_gen3,feature);

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
     accuracy(ii,1) = sum(correct)/total;
     %error(ii,1) = (1 - accuracy);
    % plot(label,'g'),hold on
    % plot(hidden)
end


Average_accuracy = sum(accuracy)/length(accuracy)


close all,
figure(5)
p1 = plot(sortgen1,1,'ro');hold on
p2 =plot(sortgen2,0,'bo');
p3 =plot(sortgen3,2,'go');

p4 = plot([threshold1 threshold1],[0 3],'y');
p5 = plot([threshold2 threshold2],[0 3],'k');
title('LDA: Three Arists different genres');
grid on
xlabel('Projected axis');
ylabel('speration of groups');
set(gca,'Fontsize',[12])
saveas(gcf,'Artists Test 1 Sample data grouping and thresholds.png')

%%

clear figure
figure(6) = figure('Position', [200, 100, 700, 500]);


subplot(1,3,1)
histogram(sortgen1), hold on
title('Projections of test 1 Artist 1');
xlabel('Smallest-to-largest projection values');
plot([threshold1 threshold1],[0 30],'r','Linewidth',3);


subplot(1,3,2)
histogram(sortgen2), hold on
title('Projections of test 1 Artist 2');
xlabel('Smallest-to-largest projection values');
plot([threshold1 threshold1],[0 30],'r','Linewidth',3);

%plot([threshold2 threshold2],[0 3],'k');

subplot(1,3,3)
histogram(sortgen3), hold on
title('Projections of test 1 Artist 3');
xlabel('Smallest-to-largest projection values');
plot([threshold2 threshold2],[0 30],'r','Linewidth',3);
saveas(gcf,'Artists Test 1 Histogram of projected values.png')


