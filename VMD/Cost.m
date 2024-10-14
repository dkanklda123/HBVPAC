function ff = Cost(X,c)
alpha = round(c(1));       % moderate bandwidth constraint：适度的带宽约束/惩罚因子
tau = 0;          % noise-tolerance (no strict fidelity enforcement)：噪声容限（没有严格的保真度执行）
K = round(c(2));              % modes：分解的模态数
DC = 0;             % no DC part imposed：无直流部分
init = 1;           % initialize omegas uniformly  ：omegas的均匀初始化
tol = 1e-7;     
%--------------- Run actual VMD code:数据进行vmd分解---------------------------
[u, u_hat, omega] = VMD(X, alpha, tau, K, DC, init, tol);
for i = 1:K
	xx= abs(hilbert(u(i,:))); %首先对分解得到的IMF分量进行希尔伯特变换，并求取幅值
	xxx = xx/sum(xx); %
    ssum=0;
    for ii = 1:size(xxx,2)
		bb = xxx(1,ii)*log(xxx(1,ii));  %最小包络熵的计算公式
        ssum=ssum+bb;  %求和运算
    end
    fitness(i,:) = -ssum;   %记着加负号！ 也是公式的一部分
end
ff = min(fitness);
end