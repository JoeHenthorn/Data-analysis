function [result,w,U,S,V,threshold1,threshold2,sortgen1,sortgen2,sortgen3] = Music_genre_trainer(train_genre1,train_genre2,train_genre3,feature)
%% The SVD
genre_1 = train_genre1;
genre_2 = train_genre2;
genre_3 = train_genre3;

A = [genre_1 genre_2 genre_3];

ngen1 = length(genre_1(1,:));
ngen2 = length(genre_2(1,:));
ngen3 = length(genre_3(1,:));
tic
[U,S,V] = svd((A),'econ');


genre_proj = S*V'; % projection onto principal components
U = U(1:feature,:);

% singularValues = diag(S);
% energy = (singularValues).^2/(sum(singularValues.^2));
% for j = 1:length(singularValues)
%     %energy_r_cum = (singularValues(1:j)).^2/(sum(singularValues.^2));
%     energy_r_cum(j) = norm(U(:,1:j)*S(1:j,1:j)*V(:,1:j)','fro').^2/norm(A,'fro').^2;
% end

% plot1 = figure(1);
% subplot(2,1,1)
% semilogy(singularValues,'b.','Markersize',[30])
% grid on 
% ylabel('\sigma (log scale)')
% %set(gca,'Fontsize',16,'Xtick',0:5:25)
% title('Singular Values for experiment ');
% 
% subplot(2,1,2)
% semilogy(energy_r_cum,'g.','Markersize',[30])
% grid on 
% title('Energies for experiment ');
% ylabel('energy (log scale)')
% set(gca,'Fontsize',12)
% %saveas(gcf,string(k) +'Singular values and energy plot semilog.png')


gen1 = genre_proj(1:feature,1:ngen1);
gen2 = genre_proj(1:feature,ngen1+1:ngen1+ngen2);
gen3 = genre_proj(1:feature,(ngen1+ngen2)+1:ngen1+ngen2+ngen3);

mean_gen1 = mean(gen1,2);
mean_gen2 = mean(gen2,2);
mean_gen3 = mean(gen3,2);

% The code below produces the scatter matrix for the Linear Discriminant Analysis

S_within = 0;

for i = 1:ngen1
    S_within = S_within + (gen1(:,i)-mean_gen1)*(gen1(:,i)-mean_gen1)';
end

for i = 1:ngen2
    S_within = S_within + (gen2(:,i)-mean_gen2)*(gen2(:,i)-mean_gen2)';
end

for i = 1:ngen3
    S_within = S_within + (gen3(:,i)-mean_gen3)*(gen3(:,i)-mean_gen3)';
end


S_between = 0;

genre_all = gen1 + gen2 + gen3;
mean_all = mean(genre_all,2);

for i = 1:ngen1
    S_between = S_between + (gen1(:,i) - mean_all)*(gen1(:,i) - mean_all)';
end

for i = 1:ngen2
    S_between = S_between + (gen2(:,i) - mean_all)*(gen2(:,i) - mean_all)';
end

for i = 1:ngen3
    S_between = S_between + (gen3(:,i) - mean_all)*(gen3(:,i) - mean_all)';
end



[V2,D] = eig(S_between,S_within);
[lambda,ind] = max(abs(diag(D)));
w = V2(:,ind); % Best lines to seperate the data
w = w/norm(w,2);

vgen1 = w'*gen1;
vgen2 = w'*gen2;
vgen3 = w'*gen3;
result = [vgen1,vgen2,vgen3];

if mean(vgen2)>mean(vgen3)
    %w = -w;
    vgen2 = -vgen2;
    flag = 1;
    vgen3 = -vgen3;
else
    flag = 0;
end
if mean(vgen1)>mean(vgen2)
    %w = -w;
    vgen1 = -vgen1;
    flag2 = 1;
    vgen2 = -vgen2;
else
    flag2 = 0;
end
% if flag == 1 ||flag2 == 1
%     w = -w;
% end


sortgen1 = -(sort(vgen1));
sortgen2 = (sort(vgen2));
sortgen3 = -(sort(vgen3));

t1 = length(sortgen1);
t2 = 1;

while sortgen1(t1)>sortgen2(t2)
    t1 = (t1 - 1);
    t2 = t2+1;
%     plot(sortgen1,0,'ro');hold on
%     plot(sortgen2,1,'bo');hold on
%     p4 = plot([threshold1 threshold1],[0 3],'k');
%     draw now 
%     if t1 == 1
%         break
%     end
end




threshold1 = (sortgen1(t1)+sortgen2(t2))/2;

t2 = length(sortgen2);
t3 = 1;

while sortgen2(t2)<sortgen3(t3)
    t2 = (t2 - 1)
    t3 = t3+1
    if t2 == 1
        break
    end
end

 threshold2 = (sortgen2(t2)+sortgen3(t3))/2;
end



