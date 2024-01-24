%function [bestFitness,record,dsp] = NRO(fobj,dim,lb,ub,maxFEs,N,F,NS)
function [bestFitness,record,dsp] = NRO(problem,dim,lb,ub,maxFEs,N,F,NS,IS,RS)
%输入:
%fobj:目标函数
%dim:维度
%lb:决策空间的下界
%ub:决策空间的上界
%maxFEs:最大函数评估次数
%N:初始种群数量
%F:步长
%NS:邻居个数
%IS:初始化策略
%RS:替换策略
%输出:
%bestFitness:目标函数最优值
%record:记录最优值历史,即[当前已使用的函数评估次数,当前所找到的目标函数最优值;...]
%X = (ub-lb)*rand(N,dim)+lb;
switch IS
    case 1
        X = (ub-lb).*rand(N,dim)+lb;%随机初始化种群
    case 2
        LHSsamples= lhsdesign(N, dim);%nomorlized training data       
        X= repmat(lb,N,1)+(repmat(ub,N,1)-repmat(lb,N,1)).*LHSsamples; %拉丁超立方随机初始化种群
    case 3
        OLHSsamples= srgtsDOEOLHS(N, dim, 'ESEA');
        X= repmat(lb,N,1)+(repmat(ub,N,1)-repmat(lb,N,1)).*OLHSsamples;%优化的拉丁超立方随机初始化种群
end
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
while maxFEs > FEs 
	for j = 1 : N
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
        switch RS
            case 1
                [~,rIdx] = max(Fitness);%替换最差的点
            case 2
                rIdx = randi(N,1);%替换随机的点
            case 3
                [~,rIdx] = min(Fitness);%替换最好的点
            case 4
                rIdx=0;%不进行替换
        end
		%[gWorstFitness,gWorstIdx] = max(Fitness);%找到最差的点
        %if UFitness <= gWorstFitness % global replacement
        %    X(gWorstIdx,:) = U; 
        %    Fitness(gWorstIdx) = UFitness; 
        %end
        if rIdx>0&&UFitness <=Fitness(rIdx)  
            X(rIdx,:)=U;
            Fitness(rIdx)=UFitness;
        end
        % record
        bestFitness = min(bestFitness,UFitness);
    end
    record = [record;[FEs,bestFitness]];
end
end