function allmode=eemd(Y,Nstd,NE)
% This is an EMD/EEMD program
%
% INPUT:
% Y: Inputted data;1-d data only
% Nstd: ratio of the standard deviation of the added noise and that of
% Y; Nstd = (0.1 ~ 0.4)*std(Y). 
% NE: Ensemble number for the EEMD, NE = 10-50.
% OUTPUT:
% A matrix of N*(m+1) matrix, where N is the length of the input
% data Y, and m=fix(log2(N))-1. Column 1 is the original data, columns 2, 3, ...
% m are the IMFs from high to low frequency, and comlumn (m+1) is the
% residual (over all trend).
%
% NOTE:
% It should be noted that when Nstd is set to zero and NE is set to 1, the
% program degenerates to a EMD program.(for EMD Nstd=0,NE=1)
% This code limited sift number=10 ,the stoppage criteria can't change.

% References:
% Wu, Z., and N. E Huang (2008),
% Ensemble Empirical Mode Decomposition: a noise-assisted data analysis method.
% Advances in Adaptive Data Analysis. Vol.1, No.1. 1-41.
%
% code writer: Zhaohua Wu.
% footnote:S.C.Su 2009/03/04
%
% There are three loops in this code coupled together.
% 1.read data, find out standard deviation ,devide all data by std
% 2.evaluate TNM as total IMF number--eq1.
% TNM2=TNM+2,original data and residual included in TNM2
% assign 0 to TNM2 matrix
% 3.Do EEMD NE times-----------loop EEMD start
% 4.add noise
% 5.give initial values before sift
% 6.start to find an IMF------IMF loop start
% 7.sift 10 times to get IMF------sift loop start and end
% 8.after 10 times sift --we got IMF
% 9.subtract IMF from data ,and let the residual to find next IMF by loop
% 6.after having all the IMFs-------------IMF loop end
% 9.after TNM IMFs ,the residual xend is over all trend
% 3.Sum up NE decomposition result--------loop EEMD end
% 10.Devide EEMD summation by NE,std be multiply back to data

%% Association: no
% this function ususally used for doing 1-D EEMD with fixed
% stoppage criteria independently.
%
% Concerned function: extrema.m
% above mentioned m file must be put together

%function allmode=eemd(Y,Nstd,NE)

%part1.read data, find out standard deviation ,devide all data by std
xsize=length(Y);
dd=1:1:xsize; 
Ystd=std(Y);
Y=Y/Ystd;    

%part2.evaluate TNM as total IMF number,ssign 0 to N*TNM2 matrix
TNM=fix(log2(xsize))-1;   % TNM=m
TNM2=TNM+2;               
for kk=1:1:TNM2,
    for ii=1:1:xsize,
        allmode(ii,kk)=0.0;  
    end
end

%part3 Do EEMD -----EEMD loop start
for iii=1:1:NE, %EEMD loop NE times EMD sum together

%part4 --Add noise to original data,we have X1
    for i=1:xsize,
        temp=randn(1,1)*Nstd; % add a random noise to Y
        X1(i)=Y(i)+temp;
    end

%part4 --assign original data in the first column
    for jj=1:1:xsize,
        mode(jj,1) = Y(jj); % assign Y to column 1of mode
    end

%part5--give initial 0to xorigin and xend
    xorigin = X1;   % 
    xend = xorigin; %

%part6--start to find an IMF-----IMF loop start
    nmode = 1;

    while nmode <= TNM,
    xstart = xend; %last loop value assign to new iteration loop
                   %xstart -loop start data
    iter = 1;      %loop index initial value

%part7--sift 10 times to get IMF---sift loop start
        while iter<=10,  
            [spmax, spmin, flag]=extrema(xstart); %call function extrema
                                                  %the usage of spline ,please see part11.
            upper= spline(spmax(:,1),spmax(:,2),dd); %upper spline bound of this sift
            lower= spline(spmin(:,1),spmin(:,2),dd); %lower spline bound of this sift
            mean_ul = (upper + lower)/2;            %spline mean of upper and lower
            xstart = xstart - mean_ul;              %extract spline mean from Xstart
            iter = iter +1;
        end

%part8--subtract IMF from data ,then let the residual xend to start to find next IMF
    xend = xend - xstart;      
    nmode=nmode+1;              

%part9--after sift 10 times,that xstart is this time IMF
    for jj=1:1:xsize,
        mode(jj,nmode) = xstart(jj);  
    end
end

%part10--after gotten all(TNM) IMFs ,the residual xend is over all trend
% put them in the last column
    for jj=1:1:xsize,
        mode(jj,nmode+1)=xend(jj);
    end

%after part 10 ,original + TNM IMFs+overall trend ---those are all in mode
    allmode=allmode+mode;   
end  %part3 Do EEMD -----EEMD loop end

%part11--devide EEMD summation by NE,std be multiply back to data
allmode=allmode/NE;    
allmode=allmode*Ystd;  
end
