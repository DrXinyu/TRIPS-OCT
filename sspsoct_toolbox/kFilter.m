function y = kFilter(x,filter)

xfft = fft(x,4096);
yfft = filter.hfft.*xfft;
yfromfft = ifft(yfft,4096);
shift = round(length(filter.states)/2);
y = circshift(yfromfft,-shift,1);

% [~,schannel] = size(x);
% 
% s = (zeros(shift,schannel));
% x = ([s;x;s]);
% 
% y = filter(Hd,x);
% y = circshift(y,-shift,1);
% y = y(shift:end-shift-1,:);











