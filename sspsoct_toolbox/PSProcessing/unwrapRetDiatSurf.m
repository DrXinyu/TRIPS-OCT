function [retuw,diatuw] = unwrapRetDiatSurf(ret,diat,init)
% takes the ret and diat of the surface signal (or any other 1D signal of
% ret/diat signal) and 'unwraps' it, starting at pixel init (defaults to
% center)

if nargin<3
    init = round(size(ret,2)/2);
end

% complex valued ret/diat vector
q = ret + 1i*diat;

% find indices that have to be wrapped
inds = mod(cumsum(cat(2,0,abs(sqrt(sum(diff(q,[],2).^2,1)))>pi)),2)>0;

if inds(init)>0
    inds = ~inds;
end

q(:,inds) = q(:,inds) - 2*pi*q(:,inds)./sqrt(sum(q(:,inds).^2,1));

retuw = real(q);
diatuw = imag(q);
