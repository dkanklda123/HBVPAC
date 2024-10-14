function [st,data] = NonlinearSig(ratio,noiseLev)

datalen = 10;% sec
fs = 600;
f_lo = 6;
f_hi = 50;
ratio = 5;
noiseLev = 0.1;
t = (0:datalen*fs-1)/fs;
sp = cos(2*pi*f_lo.*t);
[sa,~,~] = awave(fs,f_hi,ratio,datalen,1);
% st =4*sp+(2+sp).*sa;
st = (2+sp).*sa;

Noise = noiseLev*randn(size(st));
data = st + Noise;  
plot(sp+2)
PowerSpectrum(data,fs)