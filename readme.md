# RNRO

这是一个MATLAB函数，实现RNRO算法，这是一种用于解决昂贵优化问题的自然启发优化算法，该算法结合了进化计算和代理模型，以高效地探索解空间。

## 描述

该算法通过引入RBF全局代理模型找到全局的一个候选解加入到数据库中，然后执行NRO算法进行邻域的搜索，当且仅当函数的评估次数耗尽，算法终止并输出昂贵优化问题的最优解。

## 使用方法

```
[bestFitness,record,dsp]=RNRO(X,problem,lb,ub,maxFEs,N,F,NS,IS);
```

## 输入

* ‘X’: 初始种群
* ‘problem’: 问题的维度
* ‘lb’: 决策空间的下界
* ‘ub’: 决策空间的上界
* ‘maxFEs’: 最大函数评估次数
* ‘N’: 初始种群大小
* ‘F’: 步长
* ‘NS’: 邻居个数
* ‘IS’: 改进策略

## 输出

* ‘bestFitness’: 找到的最佳适应度值
* ‘record’: 优化过程中最佳适应度的历史记录
* ‘dsp’: 所使用的优化策略的描述

## 示例

```
%  定义优化问题和参数
problem=‘Ackley'
dim=10;
lb=-20;
ub=20;
maxFEs=1000;
N=50;
F=0.4;
NS=dim+3;
IS=2;

% 初始种群
LHSsamples= lhsdesign(N, D);  
X= repmat(lb,N,1)+(repmat(ub,N,1)-repmat(lb,N,1)).*LHSsamples;

% 运行RNRO算法
[bestFitness,record,dsp]=RNRO(X,problem,lb,ub,maxFEs,N,F,NS,IS);
```

该算法提供了一个框架，用于高效地解决昂贵优化问题。

注意：确保已安装必要的代理模型工具包（例如 srgts），以使该算法能够正确运行。

