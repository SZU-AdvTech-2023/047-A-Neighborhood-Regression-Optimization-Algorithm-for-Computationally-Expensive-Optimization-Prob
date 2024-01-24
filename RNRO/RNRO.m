function [bestFitness,record,dsp] = RNRO(X,problem,dim,lb,ub,maxFEs,N,F,NS,IS)
%输入:
%X:初始种群
%fobj:目标函数
%dim:维度
%lb:决策空间的下界
%ub:决策空间的上界
%maxFEs:最大函数评估次数
%N:初始种群数量
%F:步长
%NS:邻居个数
%输出:
%bestFitness:目标函数最优值
%record:记录最优值历史,即[当前已使用的函数评估次数,当前所找到的目标函数最优值;...]
%X = (ub-lb)*rand(N,dim)+lb;
%X = (ub-lb).*rand(N,dim)+lb;%初始化种群
%Fitness = arrayfun(@(j)fobj(X(j,:)),1:N)';
Fitness=arrayfun(@(j)expensive_benchmark_func(X(j,:),problem),1:N)';%计算目标函数真实值
FEs = N;
% find out the best solution 
bestFitness = min(Fitness);
record = [FEs,bestFitness];
dsp = {
	'NRO with global best strategy';
	'fixed step';
};
l=sqrt(10^(-6)*dim);
while maxFEs > FEs 
    switch IS  
       case 2
            options  = srgtsRBFSetOptions(X,Fitness);
            surrogate = srgtsRBFFit(options);
            newpoint=JADE([lb;ub],surrogate,100,100,dim,'srgtsRBFEvaluate');
            UFitness=expensive_benchmark_func(newpoint,problem);
            FEs=FEs+1;
            mind=min(pdist2(newpoint,X));
            if mind>l
                X=[X;newpoint];
                Fitness=[Fitness;UFitness];
                N=N+1;
            end
        case 3
            options  = srgtsPRSSetOptions(X,Fitness,4,'Full');
            surrogate = srgtsPRSFit(options);
            newpoint=JADE([lb;ub],surrogate,100,100,dim,'srgtsPRSEvaluate');
            UFitness=expensive_benchmark_func(newpoint,problem);
            FEs=FEs+1;
            mind=min(pdist2(newpoint,X));
            if mind>l
                X=[X;newpoint];
                Fitness=[Fitness;UFitness];
                N=N+1;
            end
        case 4
            newpoint=[];
            for i=1:dim
                ind=randi(N,1,5);
                datax=X(ind,i);
                datay=Fitness(ind);
                options  = srgtsPRSSetOptions(datax,datay,4,'Full');
                surrogate = srgtsPRSFit(options);
                lu=[min(datax);max(datax)];
                xi=JADE(lu,surrogate,100,100,1,'srgtsPRSEvaluate');
                newpoint=[newpoint,xi];
            end
            UFitness=expensive_benchmark_func(newpoint,problem);
            FEs=FEs+1;
            mind=min(pdist2(newpoint,X));
            if mind>l
                X=[X;newpoint];
                Fitness=[Fitness;UFitness];
                N=N+1;
            end
    end
   for j=1:N
   [~,sortedIdx] = sort(pdist2(X(j,:),X));%计算种群中所有点与当前点xj的欧氏距离
    NeighborIdx = sortedIdx(1:NS);%筛选距离最小的前NS个点
    NeighborX = X(NeighborIdx,:);
    NeighborFitness = Fitness(NeighborIdx);
    [~,baseIdx] = min(Fitness);%寻找种群中最优值的对应点
    baseX = X(baseIdx,:);
    increase = (max(NeighborX) - min(NeighborX))/2;%更新方向
    NeighborX = [NeighborX ones(NS,1)];
    coef = (NeighborX'*NeighborX)\NeighborX'*NeighborFitness;%B=(X'X)^{-1}X'y
    coef = coef(1:dim)';
    increase(coef>0) = -increase(coef>0);
    U = baseX + increase * F; %产生子代并检查是否越界
    ubx=(ub + baseX)/2;
    lbx=(lb+baseX)/2;
    U(U>ub) = ubx(U>ub);
    U(U<lb) = lbx(U<lb);
    %UFitness = fobj(U);
    UFitness=expensive_benchmark_func(U,problem);
    FEs = FEs + 1;
    % replacement 
    [gWorstFitness,gWorstIdx] = max(Fitness);%找到最差的点
    if UFitness <= gWorstFitness % global replacement
       X(gWorstIdx,:) = U; 
       Fitness(gWorstIdx) = UFitness; 
    end
    % rIdx =randi(N,1);%找到随机的点
    % RFitness=Fitness(rIdx);
    % if UFitness <= RFitness % global replacement
    %    X(rIdx,:) = U; 
    %    Fitness(rIdx) = UFitness; 
    % end
    % record
    bestFitness = min(bestFitness,UFitness);
   end
   %F=0.4-0.2*(FEs/maxFEs);
    record = [record;[FEs,bestFitness]];
end