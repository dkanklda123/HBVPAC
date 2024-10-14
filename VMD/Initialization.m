function [particle, GlobalBest,varargout] =  Initialization(Params,CostFunction,name)
 
nPop = Params.nPop;
VarMin = Params.VarMin;
VarMax = Params.VarMax;
VarSize = Params.VarSize;
 
%% Initialization
switch name
 
     
    % 麻雀个体    
    case 'SSA'
        
        % 捕食者个体占比
        PredatorRate = 0.5;
        % 警觉者占比
        SDRate = 0.3;
        empty_particle.Position=[];
        empty_particle.Cost=[];
 
        % 捕食者和加入者
        PredatorNumber = floor(nPop * PredatorRate);
        particle=repmat(empty_particle,nPop ,1);
 
        % 警觉者
        SDNumber = floor(nPop * SDRate);
        SD = repmat(empty_particle,SDNumber,1);
 
 
        GlobalBest.Cost=inf;
        GlobalWorst.Cost = -inf;
 
 
        % 初始化
        disp(['***初始化***'])
        for i = 1:nPop
            particle(i).Position = unifrnd(VarMin,VarMax,VarSize);
            particle(i).Cost = CostFunction(particle(i).Position);
            if GlobalBest.Cost > particle(i).Cost
                GlobalBest = particle(i);
            end
            if GlobalWorst.Cost < particle(i).Cost
                GlobalWorst = particle(i);
            end
        end
 
        % 警觉者初始化
        for i = 1:SDNumber
            SD(i).Position = unifrnd(VarMin,VarMax,VarSize);
            SD(i).Cost = CostFunction(SD(i).Position);   
        end
 
        % 挑选捕食者和加入者
        [~,index] = sort([particle.Cost]);
 
        Predator = particle(index(1:PredatorNumber));
        Joiner = particle(index(PredatorNumber+1:end));
        format long;
        disp(['初始化后的最小包络熵为：',num2str(GlobalBest.Cost),'，初始化后的最佳参数为：[',num2str(round(GlobalBest.Position)),']'])
        
end 
    
 
%%  输出
switch name
    case 'SSA'
        varargout{1} = SD;
        varargout{2} = GlobalWorst;
        varargout{3} = Predator;
        varargout{4} = Joiner;
    otherwise
       % varargout{1:4} = []; 
end
end