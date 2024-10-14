clc;clear;close all;
%% initialize parameters
fs=600;
noshufs = 100;
threshold = 0.05;
min_cycles=10;
nobins=20;
x_lims = [0 0.4]; x_bins = 20;   
y_lims = [0 6]; y_bins = 30; % for simulation data
spacing = 'linear';units='Hz';
FileName = 'Simulation_data';

%% load RR intervals time series
RR = load('RR.mat').RR; % directly load existing RR intervals

% or use synthetic ECG (following code) to aquire RR intervals
% and further processing is needed (Refer to 'Method' section in our manuscript)

% N = 200;
% Anoise = 0;
% hrmean = 75;
% hrstd = 10;
% lfh =  1;
% sfint = fs;
% [ecg, ipeaks] = ecgsyn(fs,N,Anoise,hrmean,hrstd,lfh,sfint);
% RR = diff(find(ipeaks'==3),1)/fs;
% RRtime = find(ipeaks'==3)/fs;
% RRtime = RRtime(1:end-1);

%% nonlinear coupled signal generation
EEG = load('NonlinearEEG3Hz.mat').EEG; % load existing data or generate by 
                                       % yourself (uncomment following code)

% % generate nonlinear synthetic signal
% datalen = 50;
% ratio = 9; % degree of nonlinearity
% fc = 3; % center frequency of the simulated signal
% [sa,~,~] = awave(fs,fc,ratio,datalen,1);
% Xa1 = sa(1:fs*datalen);
% 
% % generate nonlinear coupled signal
% Sp = (RR-median(RR))/max(RR-median(RR));
% st = (0.75*(1+Sp)+0.25).*Xa1;
% noiseLev = 0.5;
% Noise = noiseLev*randn(size(st));
% EEG = st + Noise;

%% Decompose RR intervals into LF HF components using VMD 
alpha = 50000;
tau = 0;
K = 3;
DC = 0;
init = 1;
tol = 1e-7;
[RRimf,~,RRomega] = VMD(RR,alpha,tau,K,DC,init,tol);
[~, sortIndex] = sort(RRomega(end,:),'descend');
RRimf = RRimf(sortIndex,:);

%% Decompose nonlinear coupled signal by VMD 
% alpha = 2000;
% tau = 0;
% K = 2;
% DC = 0;
% init = 1;
% tol = 1e-7;
% [EEGimf,~,EEGomega] = VMD(EEG,alpha,tau,K,DC,init,tol);
% [~, sortIndex] = sort(EEGomega(end,:),'descend');
% EEGimf = EEGimf(sortIndex,:);

%% Decompose nonlinear coupled signal by EEMD
% oth = 0.3;
% EEGimf = feemd(EEG, 0.1, 2, -1); 
% [datalen,modenum] = size(EEGimf);
% if datalen<modenum
%     EEGimf = EEGimf';
% end
% EEGimf(:,end)=[]; % take off trend
% % orthogonal check
% o1=1;
% while o1<size(EEGimf,2)
%     or_up=sum(EEGimf(:,o1).*EEGimf(:,o1+1));
%     or_down1=sqrt(sum(EEGimf(:,o1).*EEGimf(:,o1)));
%     or_down2=sqrt(sum(EEGimf(:,o1+1).*EEGimf(:,o1+1)));
%     or_idx=or_up/or_down1/or_down2;
%     if abs(or_idx)>oth %
%         EEGimf(:,o1)=EEGimf(:,o1)+EEGimf(:,o1+1);
%         EEGimf(:,o1+1:end-1)=EEGimf(:,o1+2:end);
%         EEGimf(:,end)=[];
%     else
%         o1=o1+1;
%     end
% end
% EEGimf = EEGimf';
%% Decompose nonlinear coupled signal by EMD
% [EEGimf, residual, info] = emd(EEG);
% [datalen,modenum] = size(EEGimf);
% if datalen<modenum
%     EEGimf = EEGimf';
% end
% EEGimf(:,end)=[]; % take off trend
% % orthogonal check
% oth=0.3;
% o1=1;
% while o1<size(EEGimf,2)
%     or_up=sum(EEGimf(:,o1).*EEGimf(:,o1+1));
%     or_down1=sqrt(sum(EEGimf(:,o1).*EEGimf(:,o1)));
%     or_down2=sqrt(sum(EEGimf(:,o1+1).*EEGimf(:,o1+1)));
%     or_idx=or_up/or_down1/or_down2;
%     if abs(or_idx)>oth %
%         EEGimf(:,o1)=EEGimf(:,o1)+EEGimf(:,o1+1);
%         EEGimf(:,o1+1:end-1)=EEGimf(:,o1+2:end);
%         EEGimf(:,end)=[];
%     else
%         o1=o1+1;
%     end
% end
% EEGimf = EEGimf';

%% 
imfs1 = [EEGimf;RRimf];
[datalen,modenum] = size(imfs1);
if datalen<modenum
    imfs1 = imfs1';
end
%% EEG vaild mode
[~,~,P,F,~,~,~] = makeHAPF(EEGimf',fs,min_cycles); % P: Phase F: Cycle-frequency
[~,EEGnum] = size(P);
EEGnum = EEGnum+1;

%% HAPF of imfs1
[~,~,P,F,~,~,~] = makeHAPF(imfs1,fs,min_cycles);
[~,modenum] = size(P);
% get instataneous frequency
instF = fs*diff(P)./(2*pi);
instF = [instF;instF(end,:)];

%% show 1st IMFs components,ratio of cycle frequency for each pairs of IMFs.
[Fratio,Fpercent] = getIMFclclen7(imfs1,instF,F);

%% get IMFs information for each pairs of IMFs
% i:phase. j:amplitude.
for i=EEGnum:modenum % get phase from amplitude
    for j=1:i-1 % amplitude
        %% available points according to cycleF ratio(Niquist rate)
        if Fpercent(j,i)>10
            pts = find(Fratio{j,i}>=2);%
        else
            pts = [];
        end
        PTS{j,i} = pts;
        if ~isempty(PTS{j,i})
            %% information of IMFs (Phase)
            [~,~,P_tempP{j,i},~,cycle_bounds_tempP{j,i},bands_tempP{j,i}]...
                = makeHAPF1(imfs1(:,i),fs,min_cycles,PTS{j,i});  %%% AAC PPC 3P 2A
            %% information of IMFs (Amplitude)
            [~,A_tempA{j,i},~,~,cycle_bounds_tempA{j,i},bands_tempA{j,i}]...
                = makeHAPF1(imfs1(:,j),fs,min_cycles,PTS{j,i});
        end
    end
end
%% Computing amplitude vs. phase distributions.
MI = zeros(modenum,modenum);
for i=EEGnum:modenum % phase
    for j=1:i-1 % amplitude
        if ~isempty(PTS{j,i}) && size(P_tempP{j,i},2)~=0 && ~isempty(A_tempA{j,i}) && size(bands_tempA{j,i},1)~=0
            %% plot and computing amplitude vs. phase distributions
            [M(:,j,i),validMode{j,i}] = PACdistribution2(P_tempP{j,i},A_tempA{j,i},nobins);
            %% Computing modulation index(with consideration of valid mode)
            if sum(validMode{j,i})~=0
                MI(j,i)=inv_entropy(M(:,j,i));
            else
                MI(j,i)=0;
            end
            % clear M validMode
        end
    end
end
%% inverse entropy p-values.
MI_z_scores = zeros(size(MI));
MI_z_thresh = zeros(size(MI));
MImean = zeros(size(MI));
MIstd = zeros(size(MI));
if noshufs>0
    for i=EEGnum:modenum % phase
        for j=1:i-1 % amplitude
            if ~isempty(PTS{j,i}) && ~isempty(A_tempA{j,i}) && size(P_tempP{j,i},2)~=0
                [MI_z_scores(j,i),MImean(j,i),MIstd(j,i)] = ...
                    inv_entropy_distn(noshufs,nobins,...
                    A_tempA{j,i},P_tempP{j,i},...
                    cycle_bounds_tempA{j,i},cycle_bounds_tempP{j,i},MI(j,i));
            end
        end
    end
    MI_z_thresh = inv_entropy_thresh(triu(MI,1),triu(MI_z_scores,1),threshold);
end
%% scatter plot and comodulogram
if max(max(MI_z_thresh))>0
    Binned_zMI = emd_inv_entropy_plot_bMI1_loop(MI_z_thresh,F,PTS,x_lims,x_bins,y_lims,y_bins,spacing,units,['Zscore_' FileName]);
else
    Binned_zMI = zeros(size(Binned_MI));
end

%% creat xlabel & ylabel
x_bin = linspace(x_lims(1),x_lims(2),x_bins+1);
y_bin = linspace(y_lims(1),y_lims(2),y_bins+1);
for m=1:x_bins+1
    xlabels{m}=num2str(x_bin(m));
end
for n=1:y_bins+1
    ylabels{n}=num2str(y_bin(n));
end

tick1 = 1:(x_bins+1); tick1 = tick1([1 9 21]);
tick2 = 1:(y_bins+1); tick2 = tick2([1 16 31]);
xlabels = xlabels([1 9 21]);
ylabels = ylabels([1 16 31]);

%% plot & save MI and Binned_MI
figure,
colorplot(imgaussfilt(Binned_zMI,2))
colormap(bone)
axis xy
set(gca,'XTick',tick1,'YTick',tick2,'XTickLabel',xlabels,'YTickLabel',ylabels);
set(gca,'FontSize',14);
title('PAC (zscore)','FontSize',14)
xlabel('RR Phase (Hz)','FontSize',14);
ylabel('EEG Amp. (Hz)','FontSize',14);



