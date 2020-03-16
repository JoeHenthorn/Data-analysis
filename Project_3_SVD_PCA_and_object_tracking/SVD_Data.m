function [singularValues, energy, energy_r_cum,plot1,plot2,plot3,plot4] = SVD_Data(Position_Data,k)


% The code below makes the vectors all start at the same place.
%Peaks_1 = peaks(Posy_1) 

Posx_1 = cell2mat(Position_Data(1,:));
Posy_1 = cell2mat(Position_Data(2,:));
Posx_2 = cell2mat(Position_Data(3,:));
Posy_2 = cell2mat(Position_Data(4,:));
Posx_3 = cell2mat(Position_Data(5,:));
Posy_3 = cell2mat(Position_Data(6,:));


[Max_1, indx] =    max(Posy_1(1:50));
t_max_1 = ind2sub(size(Posy_1),indx);
Posx_1 = Posx_1(t_max_1:end);
Posy_1 = Posy_1(t_max_1:end);

[Max_2, indx] =    max(Posy_2(1:50));
t_max_2 = ind2sub(size(Posy_2),indx);
Posx_2 = Posx_2(t_max_2:end);
Posy_2 = Posy_2(t_max_2:end);

[Max_3, indx] =    max(Posy_3(1:50));
t_max_3 = ind2sub(size(Posy_3),indx);
Posx_3 = Posx_3(t_max_3:end);
Posy_3 = Posy_3(t_max_3:end);

%%% The code below makes the vectors all have same length by truncating the
%%% ends of the two longer vector with the shortest vector.
L1 = length(Posx_1);
L2 = length(Posx_2);
L3 = length(Posx_3);
L_all = [L1 L2 L3];
trunc = min(L_all);

Posx_1 = Posx_1(1:trunc);
Posx_2 = Posx_2(1:trunc);
Posx_3 = Posx_3(1:trunc);

Posy_1 = Posy_1(1:trunc);
Posy_2 = Posy_2(1:trunc);
Posy_3 = Posy_3(1:trunc);

L1 = length(Posy_1);
L2 = length(Posy_2);
L3 = length(Posy_3);
L_all_y = [L1 L2 L3];

%%% Plot of camera recording vert_dist
plot1 = figure(1), hold on, set(gcf, 'Position', [200, 100, 700, 500]);

plot(1:length(Posy_1),Posy_1,'g','Linewidth',2)
title("Spring motion from three cameras. Experiment #" + string(k))
xlabel('Time');
ylabel('Vertical displacement');
set(gca,'Fontsize',12)
xlim([1 length(Posy_1)]);
grid on 

plot(1:length(Posy_2),Posy_2,'Linewidth',2)
xlabel('Time');
ylabel('Vertical displacement');
grid on 
xlim([1 length(Posy_2)]);
set(gca,'Fontsize',12)

plot(1:length(Posy_3),Posy_3,'r','Linewidth',2);
xlabel('Time');
ylabel('Vertical displacement');
grid on 
xlim([1 length(Posy_3)]);
set(gca,'Fontsize',12)
legend('Camera 1','Camera 2','Camera 3')
saveas(gcf,string(k) +'Plot of camera recording vert_dist.png')

% POD / PCS / SVD 

A = [Posx_1;Posy_1;Posx_2;Posy_2;Posx_3;Posy_3];

[U,S,V] = svd(A,'econ');

singularValues = diag(S);
energy = (singularValues).^2/(sum(singularValues.^2));
for j = 1:length(singularValues)
    %energy_r_cum = (singularValues(1:j)).^2/(sum(singularValues.^2));
    energy_r_cum(j) = norm(U(:,1:j)*S(1:j,1:j)*V(:,1:j)','fro').^2/norm(A,'fro').^2;

end
plot2 = figure(2);
subplot(2,1,1)
semilogy(singularValues,'b.','Markersize',[30])
grid on 
ylabel('\sigma (log scale)')
%set(gca,'Fontsize',16,'Xtick',0:5:25)
title('Singular Values for experiment ' + string(k));

subplot(2,1,2)
semilogy(energy_r_cum,'g.','Markersize',[30])
grid on 
title('Energies for experiment '+ string(k));
ylabel('energy (log scale)')
set(gca,'Fontsize',12)
saveas(gcf,string(k) +'Singular values and energy plot semilog.png')

% Data from each channel & SVD approximation ranked 1 to 6
figure(3), set(gcf, 'Position', [100, 100, 1050, 1150]); 
for j=1:6
    
     subplot(6,2,(2*j)-1)
%     if j == 1
%         col1 = 'g';
%     elseif j == 3
%         col2 = 'r';
%     elseif j == 5
%         col3 = 'k';
%     end
%     if j == 1 || j == 2
%         c = col1;
%     elseif j == 3 || j == 4
%         c = col2;
%     elseif j == 5 || j == 6
%         c = col3;
    c = {'g' ;'g' ; 'r' ;'r' ; 'k';'k'};
    c = char(c(j));
    plot(A(j,:),c,'Linewidth',2)
    
    if j == 1
        title('Data X_a. Ex: '+ string(k));
        ylabel('Horizontal disp.')
    elseif j == 2
        title('Data Y_a. Ex: '+ string(k));
        ylabel('Vertical disp.')
    elseif j == 3
        title('Data X_b. Ex: '+ string(k));
        ylabel('Horizontal disp.')
    elseif j == 4
        title('Data Y_b. Ex: '+ string(k));
        ylabel('Vertical disp.')
    elseif j == 5
        title('Data X_c. Ex: '+ string(k));
        ylabel('Horizontal disp.')
    elseif j == 6
        title('Data Y_c. Ex: '+ string(k));
        ylabel('Vertical disp.')
    end
    
    xlabel('time')
    xlim([1 length(A)]);
    grid on 
    set(gca,'Fontsize',10)
end

system_1 = [];
for  j=1:6
    system1 = U(:,j)*S(j,j)*V(:,j)';
    system_1 = [system_1; sum(system1)];
    subplot(6,2,2*j)
    
    plot(system1(j,:),'Linewidth',2)
    if j == 1
        title('Rank 1-6 PCA individual approximations. Ex '+ string(k));
    end
    ylabel('Vertical disp.')
    xlabel('time')
    xlim([1 length(system1)]);
    set(gca,'Fontsize',10)
    grid on 
    
end
saveas(gcf,string(k) +'Data from each channel & SVD approximation ranked 1 to 6.png')

%% Average data vs rank 1, 2, 3 approximations
plot3 = figure(4), set(gcf, 'Position', [200, 100, 700, 700]);  
subplot(4,1,1)
A_avg = sum(A)/length(A(1,:));
plot(A_avg,'r','Linewidth',4)
grid on 
title('Average position of can for experiment '+ string(k))
ylabel('Vertical displacement')
xlabel('time')
xlim([1 length(A_avg)]);
set(gca,'Fontsize',12)

subplot(4,1,2)
rank_1 = system_1(1,:)/length(system1(1,:));
plot(rank_1,'b','Linewidth',4)
title('PCA: Rank 1 approximation for experiment '+ string(k))
ylabel('Vertical displacement')
xlabel('time')
xlim([1 length(rank_1)]);
grid on 
set(gca,'Fontsize',12)

subplot(4,1,3)
rank_2 = sum(system_1(1:2,:))/length(system1(1,:));
plot(rank_2,'b','Linewidth',4)
title('PCA: Rank 1-2 Avg approximation for experiment '+ string(k))
ylabel('Vertical displacement')
xlabel('time')
xlim([1 length(rank_2)]);
grid on 
set(gca,'Fontsize',12)

subplot(4,1,4)
rank_3 = sum(system_1(1:3,:))/length(system1(1,:));
plot(rank_3,'b','Linewidth',4)
title('PCA: Rank 1-3 avg approximation for experiment '+ string(k))
ylabel('Vertical displacement')
xlabel('time')
xlim([1 length(rank_3)]);
grid on 
set(gca,'Fontsize',12)
saveas(gcf,string(k) +'Average data vs rank 1, 2, 3 approximations.png')

%%
plot4 = figure(5)
proj_mat = U'*A;
plot(proj_mat(1,:),'b','Linewidth',4)


end

