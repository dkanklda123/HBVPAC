function [s,am,ph] = awave(fs,f,cycleR,datalen,option)
% This program generate asymetical (nonlinear) signal.
%-------------------------------------------------------------------------
% INPUT : 
%        fs : sampling frequency(Hz).
%        f : target frequency(Hz).
%        cycleR : cycle ratio(T1/T2).
%        datalen : data length(sec).
% OUTPUT :
%        s : asymetical synthetic signal.
%        am : instataneous amplitude.
%        ph : instataneous phase.
%-------------------------------------------------------------------------
% Author : Chien-Hung Yeh
%-------------------------------------------------------------------------
% decide periods of tempate
T = 1/f; % target period
T1 = T*cycleR/(cycleR+1)*2; % larger period 
T2 = T*1/(cycleR+1)*2; % smaller period

% make up template
t=[0:T1*fs/4-1]/fs; % 1/4 of T1
p1=cos(2*pi*t/T1-pi/2); % phase of cosine:[-pi/2 0]
t=[0:T2*fs/2-1]/fs; % 1/2 of T2
p2=cos(2*pi*t/T2); % phase of cosine:[0 pi]
t=[0:T1*fs/4-1]/fs; % 1/4 of T1
p3=cos(2*pi*t/T1-pi); % phase of cosine:[-pi -pi/2]
s0=[p1,p2,p3];

% repeat template
n = ceil(datalen/length(s0)*fs);
s=[];
for i=1:n
    s=[s,s0]; % nonstationary 4*rand(1)*
end
    
% time & inst phase & inst amplitude 
t=[0:length(s)-1]/fs;
hs=hilbert(s);
ph=phase(hs);
ph=mod(ph,2*pi);
am=abs(hs);

% amplitude normalization
s=s./mean(am);
hs=hilbert(s);
ph=phase(hs);
ph=mod(ph,2*pi);
am=abs(hs);

if option == 0
    figure
    subplot(3,1,1);plot(t,s,'k','linewidth',2);
%     axis tight;
    set(gca,'FontSize',14);
    set(gca,'box','off','linewidth',2);
    title(['Target Frequency is ' num2str(f) 'Hz (Ratio of Period is ' num2str(cycleR) ')'],'FontSize',14)
    subplot(3,1,2);plot(t,am,'k','linewidth',2);axis tight;
    set(gca,'FontSize',14);
    set(gca,'box','off','linewidth',2);
    title('Instataneous Amplitude','FontSize',14)
    subplot(3,1,3);plot(t,ph,'k','linewidth',2);axis tight;
    set(gca,'FontSize',14);
    set(gca,'box','off','linewidth',2);
    title('Instataneous Phase (radian)','FontSize',14)
    xlabel('Time(Sec)','FontSize',14)
end

% % plot for components
% tt1_1=[0:T1*fs/4]/fs; % 1/4 of T1
% pp1_1=cos(2*pi*tt1_1/T1-pi/2); % phase of cosine:[-pi/2 0]
% tt1_2=[T1*fs/4:3*T1*fs/4]/fs; % 1/2 of T1
% pp1_2=cos(2*pi*tt1_2/T1-pi/2); % phase of cosine:[-pi/2 0]
% tt1_3=[3*T1*fs/4:T1*fs]/fs; % 1/4 of T1
% pp1_3=cos(2*pi*tt1_3/T1-pi/2); % phase of cosine:[-pi/2 0]
% figure
% subplot(311)
% plot(tt1_1-T1/2,pp1_1,'k','linewidth',4);hold on;
% plot(tt1_2-T1/2,pp1_2,'k--','linewidth',4);
% plot(tt1_3-T1/2,pp1_3,'k','linewidth',4);
% axis tight;
% set(gca,'FontSize',14);
% set(gca,'box','off','linewidth',2);
% title('Period is 0.08 Sec (4/5 of Target Period)','FontSize',14)
% xlabel('Time(Sec)','FontSize',14)
% % plot for components
% tt2_1=[0:T2*fs/4]/fs; % 1/4 of T1
% pp2_1=cos(2*pi*tt2_1/T2-pi/2); % phase of cosine:[-pi/2 0]
% tt2_2=[T2*fs/4:3*T2*fs/4]/fs; % 1/2 of T1
% pp2_2=cos(2*pi*tt2_2/T2-pi/2); % phase of cosine:[-pi/2 0]
% tt2_3=[3*T2*fs/4:T2*fs]/fs; % 1/4 of T1
% pp2_3=cos(2*pi*tt2_3/T2-pi/2); % phase of cosine:[-pi/2 0]
% subplot(312)
% plot(tt2_1-T2/2,pp2_1,'r--','linewidth',4);hold on;
% plot(tt2_2-T2/2,pp2_2,'r','linewidth',4);
% plot(tt2_3-T2/2,pp2_3,'r--','linewidth',4);
% axis tight;
% set(gca,'FontSize',14);
% set(gca,'box','off','linewidth',2);
% title('Period is 0.02 Sec (1/5 of Target Period)','FontSize',14)
% xlabel('Time(Sec)','FontSize',14)
% % plot
% subplot(313)
% plot(tt1_1-T1/4-T2/4,pp1_1,'k','linewidth',4);hold on;
% plot(tt2_2-T2/2,pp2_2,'r','linewidth',4);
% plot(tt1_3-3*T1/4+T2/4,pp1_3,'k','linewidth',4);
% axis tight;
% set(gca,'FontSize',14);
% set(gca,'box','off','linewidth',2);
% title('Period is 0.1 Sec (Target Period)','FontSize',14)
% xlabel('Time(Sec)','FontSize',14)