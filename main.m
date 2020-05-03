clc,clear,close all
warning off
format longG

c1 = 1.4995;  
c2 = 1.4995;
Vmin = -1;
Vmax = 1;
maxgen = 200;   % 迭代次数
sizepop = 20;  % 种群数量
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
beta = 1.5; % beta控制常数
nvar = 2;   % 2个未知量
% 初始化种群
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
    V(i,1) = 0;
    V(i,2) = 0;
end
% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % 全局最佳
gbest=pop;                % 个体最佳
fitnessgbest=fitness;     % 个体最佳适应度值
fitnesszbest=bestfitness; % 全局最佳适应度值

% 迭代寻优
for iter=1:maxgen
   
    for j=1:sizepop
        mbest = mean(gbest);
        % 速度更新
        r1 = rand;
        r2 = rand;
        V(j,:) = ( c1*r1*( gbest(j,:) ) + c2*r2*( zbest ) ) / ( c1*r1+c2*r2 );
        % V--x1
        if V(j,1)>Vmax
            V(j,1)=Vmax;
        end
        if V(j,1)<Vmin
            V(j,1)=Vmin;
        end
        % V--x2
        if V(j,2)>Vmax
            V(j,2)=Vmax;
        end
        if V(j,2)<Vmin
            V(j,2)=Vmin;
        end
        
        % 个体更新
        u = rand;
        if u>0.5
            pop(j,:) = Levy(beta,nvar).*(pop(j,:)-zbest) + ( 0.5+0.5*(maxgen-iter)/maxgen ).*abs(mbest-pop(j,:)) .*log(u);
%             pop(j,:) = Levy(beta,nvar) + ( 0.5+0.5*(maxgen-iter)/maxgen ).*abs(mbest-pop(j,:)) .*log(u);
        else
            pop(j,:) = Levy(beta,nvar).*(pop(j,:)-zbest) - ( 0.5+0.5*(maxgen-iter)/maxgen ).*abs(mbest-pop(j,:)) .*log(u);
%             pop(j,:) = Levy(beta,nvar) - ( 0.5+0.5*(maxgen-iter)/maxgen ).*abs(mbest-pop(j,:)) .*log(u);
        end
        
        % x1
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % 适应度更新
        fitness(j) = fun(pop(j,:));
        
        % 比较  个体间比较
        if fitness(j)<fitnessgbest(j)
            fitnessgbest(j) = fitness(j);
            gbest(j,:) = pop(j,:);
        end
        if fitness(j)<bestfitness
            bestfitness = fitness(j);
            zbest =  pop(j,:);
        end
        
    end
    fitness_iter(iter) = bestfitness;
   
end
disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
% plot(fitness_iter,'ro-','linewidth',2)
loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on