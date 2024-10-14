function [ GlobalBest,BestCost] =SSA (particle,GlobalBest,GlobalWorst,SD,Predator,Joiner,Params,CostFunction)
 
MaxIter = Params.MaxIter;
nPop = Params.nPop;
VarMin = Params.VarMin;
VarMax = Params.VarMax;
VarSize = Params.VarSize;
nVar=2;%size(VarSize,2);
BestCost = zeros(1,MaxIter);
 
%% Main loop
for i = 1:MaxIter
    for j = 1:length(Predator)
        alarm =  randn ;
        ST = randn;
        if alarm < ST
            Predator(j).Position = Predator(j).Position .* exp( -j /MaxIter);
        else
            Predator(j).Position = Predator(j).Position + randn * ones(VarSize);
        end
        Predator(j).Position = max(VarMin,Predator(j).Position);
        Predator(j).Position = min(VarMax,Predator(j).Position);
        Predator(j).Cost = CostFunction(Predator(j).Position);
      
    end
      [~,idx] = min([Predator.Cost]);
      BestPredator = Predator(idx);
    % 加入者更新
    for j = 1: nPop - length(Predator)
        if j + length(Predator)> nPop/2
            Joiner(j).Position =  randn .* exp( (GlobalWorst.Position - Joiner(j).Position) / j^2);
            
        else
            A = randi([0,1],1,nVar);
            A(~A) = -1;
            Ahat = A' / (A * A');
                Joiner(j).Position = BestPredator.Position + abs(Joiner(j).Position - BestPredator.Position) * Ahat * ones(VarSize);  
        end
        Joiner(j).Position = max(VarMin,Joiner(j).Position);
        Joiner(j).Position = min(VarMax,Joiner(j).Position);
        Joiner(j).Cost = CostFunction(Joiner(j).Position);
    end
    
    % 警觉者更新
    for j = 1:length(SD)
        if SD(j).Cost > GlobalBest.Cost
            SD(j).Position = GlobalBest.Position + randn * abs( SD(j).Position - GlobalBest.Position);
        
        elseif SD(j).Cost == GlobalBest.Cost
            SD(j).Position = SD(j).Position + (rand*2-1) * (abs( SD(j).Position - GlobalWorst.Position)./ ((SD(j).Cost - GlobalWorst.Cost) + 0.001));
        end
        SD(j).Position = max(VarMin,SD(j).Position);
        SD(j).Position = min(VarMax,SD(j).Position);
        SD(j).Cost = CostFunction(SD(j).Position);
    end
  
    
    
% 更新
particle = [Predator;Joiner;SD];
for m = 1:length(particle)
    if GlobalBest.Cost > particle(m).Cost
        GlobalBest = particle(m);
    end
    if GlobalWorst.Cost < particle(m).Cost
        GlobalWorst = particle(m);
    end
end
 
BestCost(i) = GlobalBest.Cost;
 
 
disp(['第',num2str(i),'次寻优的最小包络熵为：',num2str(BestCost(i)),'，对应最佳参数为：[',num2str(round(GlobalBest.Position)),']'])
   
 
end
end