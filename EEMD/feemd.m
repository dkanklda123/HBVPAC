function [allmode] = feemd(Y, NoiseLevel, NE, TNM, varargin)
%---------------------------------------------------------------
% INPUT : 
%        Y : input signal
%        NoiseLevel : level of added noise
%        NE : number of ensemble
%        TNM : number of required imf; if it is less than zero, then automatically determine the number
%
%-------------Additional Input Properties-----------------------------------------
%
%        toFlip : 0=> Original EEMD, References[2] ; 1=> Add anti-phase noise into signal, the way is the same as CEEMD, References[3]
%        numIteration : number of sifting iteration
%        typeSpline : 1=> clamped spline; 2=> not a knot spline;
%        toModify : 0=> None ; 1=> Apply modified linear extrapolation to boundary ; 2 => Mirror Boundary
%        randType : 1=> uniformly distributed white noise; 2=> gaussian white noise
%        seedNo : random seed used for white noise; The value of seed must be an integer between 0 and 2^32 - 1
%        checkSignal : 1=> verify if input signal had NaN or Infinity elements; 0=> Not verify input signal 
%
% OUTPUT :
%        allmode : returned imf
%
% References
%    [1] N. E. Huang, Z. Shen, S. R. Long, M. L. Wu, H. H. Shih, Q. Zheng, N. C. Yen, C. C. Tung, and H. H. Liu, 
%        "The Empirical Mode Decomposition and Hilbert spectrum for nonlinear and non-stationary time series analysis,"
%        Proc. Roy. Soc. London A, Vol. 454, pp. 903¡V995, 1998.
%    [2] Z. Wu and N. E. Huang, "Ensemble Empirical Mode Decomposition: A
%        Noise-Assisted Data Analysis Method," Advances in Adaptive Data Analysis, vol. 1, pp. 1-41, 2009.
%    [3] J. R. Yeh, J. S. Shieh and N. E. Huang, "Complementary ensemble empirical
%        mode decomposition: A novel noise enhanced data analysis method," 
%        Advances in Adaptive Data Analysis, vol. 2, Issue 2, pp. 135¡V156, 2010.
%
%  Copyright (C) RCADA, National Central University; 2013
%  ***** For Educational and Academic Purposes Only ********* 
%  Author : Yung-Hung Wang, RCADA, NCU
%


%fprintf('Copyright (C) RCADA, NCU.\n');

allmode = [];
%verifyCode = 20041017; % For verification of emd

[Y,NoiseLevel,NE,TNM,toFlip,numIteration,typeSpline,toModifyBC,randType,seedNo,IsInputOkay] = parse_checkProperty(Y, NoiseLevel, NE, TNM, varargin);

if(~IsInputOkay)
    fprintf('ERROR : The process is not executed.\n');
    return;
end

if (NoiseLevel == 0)
    allmode = rcada_emd(Y, toModifyBC, typeSpline, TNM, numIteration);
    return;
end

xsize = size(Y,2);
Ystd = std(Y);	

allmode = zeros(xsize,TNM);
imf = zeros(xsize,TNM);

savedState = set_seed(seedNo);
	
if (toFlip)
   NE = 2*NE; % YHW0202_2011: flip noise to balance the perturbed noise
end

for iii=1:NE  % ensemble loop    
    if (toFlip)
        if (mod(iii,2) ~= 0)
            if (randType == 1) % White Noise
               temp = ((2*rand(1,xsize)-1)*NoiseLevel).*Ystd; 
            elseif (randType == 2) % Gaussian Noise
               temp = (randn(1,xsize)*NoiseLevel).*Ystd; 
            end
         else % Even number
           temp = -temp;
        end
    else % toFlip = 0
        if (randType == 1)
            temp = (2*rand(1,xsize)-1)*NoiseLevel.*Ystd; % temp is Ystd*[0 1]	 
        elseif (randType == 2)
            temp = randn(1,xsize)*NoiseLevel.*Ystd; % temp is Ystd*[0 1]
        end 
    end
    xend =  Y + temp;
    imf = rcada_emd(xend, toModifyBC, typeSpline, TNM, numIteration);
    allmode = allmode + imf; 
end % iii: ensemble loop

return_seed(savedState);
allmode = allmode/NE;
return; % end eemd

end

function savedState = set_seed(seedNo)
defaultStream = RandStream.getGlobalStream;
% defaultStream = RandStream.getGlobalStream;
savedState = defaultStream.State;
rand('seed',seedNo);
randn('seed',seedNo);

end

function return_seed(savedState)
RandStream.getDefaultStream.State = savedState;
end

function [Y, NoiseLevel, NE, TNM, toFlip, numIteration, typeSpline,toModifyBC,randType,seedNo, IsInputOkay] = parse_checkProperty(Y, NoiseLevel, NE, TNM, varargin)
% Default Parameters
toFlip = 0; % Original EEMD
numIteration = 10; % numIteration = 10
typeSpline = 2;
toModifyBC = 1;
randType = 2;
seedNo = now;
checkSignal = 0;
IsInputOkay = true;

if(~isempty(varargin{1}))

for iArg = 1 : length(varargin{1});
    
if(iArg == 1)
   toFlip = varargin{1}{iArg};
   if(toFlip ~= 0 && toFlip ~= 1)
    fprintf('ERROR : toFlip must be 0 (Off) or 1 (On).\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 2)
   numIteration = varargin{1}{iArg};
   if(numIteration < 1 || (mod(numIteration, 1) ~= 0))
    fprintf('ERROR : Number of Iteration must be an integer more than 0.\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 3)
    typeSpline = varargin{1}{iArg};
    if(typeSpline ~= 1 && typeSpline ~= 2 && typeSpline ~= 3)
    fprintf('ERROR : typeSpline must be 1 (clamped spline); 2 (not a knot spline).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 4)
    toModifyBC = varargin{1}{iArg};
    if(toModifyBC ~= 0 && toModifyBC ~= 1 && toModifyBC ~= 2)
    fprintf('ERROR : toModifyBC must be 0 (None) ; 1 (modified linear extrapolation); 2 (Mirror Boundary)\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 5)
    randType = varargin{1}{iArg};
    if(randType ~= 1 && randType ~= 2)
    fprintf('ERROR : randType must be 1 (uniformly distributed white noise) ; 2 (gaussian white noise).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 6)
    seedNo = varargin{1}{iArg};
    if(seedNo < 0 || seedNo >  2^32-1 || (mod(seedNo, 1) ~= 0))
    fprintf('ERROR : The value of seed must be an integer between 0 and 2^32 - 1. \n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 7)
    checkSignal = varargin{1}{iArg};
    if(checkSignal ~= 0 && checkSignal ~= 1)
    fprintf('ERROR : Number of checksignal must be 1 (Yes) or 0 (No).\n');
    IsInputOkay = false;
    return;
    end
end

end

end

if(NoiseLevel == 0)
    fprintf('If NoiseLevel is ZERO, EEMD algorithm will be changed to EMD algorithm.\n');
end
if ((NE < 1) || (mod(NE, 1) ~= 0))
    fprintf('ERROR : Number of Ensemble must be integer more than 0.\n');
    IsInputOkay = false;
    return;
end


[m,n] = size(Y);
if(m ~= 1)
    if((n ~= 1))
       fprintf('ERROR : EMD could not input matrix array !\n');
       IsInputOkay = false;
       return;
    else
        Y =Y';
        xsize = m;
    end
else
    xsize = n;
end

if (checkSignal == 1)
    if((any(isinf(Y(:)) == 1)) || (any(isnan(Y(:)) == 1)))
        fprintf('ERROR : The input signal has NaN or Infinity elements.\n');
        IsInputOkay = false;
        return;
    end
end

if(mod(TNM, 1) ~= 0)
    fprintf('ERROR : TNM must be an integer more than 0. \n');
    IsInputOkay = false;
    return;
end

if (TNM <= 0) % automatic estimating number of imf 
    TNM=fix(log2(xsize));
end

end
