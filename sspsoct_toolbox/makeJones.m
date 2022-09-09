function Jout = makeJones(r,d)
%Jout = makeJones(r,d) generates the general Jones matrix, defined by the
%retaration 'r' and diattenuation 'd'. 'r' and 'd' both have 3 components
%in the first dimension, but can by an array of any matching dimension in
%the additional dimensions.
%If r is only a three component vector, Jout will be reshape to 2x2
%elements, otherwise Jout is 4 along the first dimension, and matches r for
%the additional dimensions.

if nargin<2 || isempty(d)
    d = zeros(size(r));
end

dim = size(r);

f = (d - 1i*r)/2;
c = sqrt(sum(f.^2));

sinch = sinh(c)./c;
sinch(c==0) = 1;
Jout = bsxfun(@times,[1;0;0;1],cosh(c(:,:))) + sinch(:,:).*(bsxfun(@times,[1;0;0;-1],f(1,:)) + bsxfun(@times,[0;1;1;0],f(2,:)) + bsxfun(@times,[0;1i;-1i;0],f(3,:)));

if numel(r)==3
    Jout = reshape(Jout,[2,2]);
else
    Jout = reshape(Jout,[4,dim(2:end)]);
end