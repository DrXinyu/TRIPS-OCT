function out = makeRot(omega,dim)
% out = makeRot(omega,dim)
% Converts omega into a linearized 3x3 rotation matrix, intepreting the
% three indices along dim as the three components of the rotation vector. 
% size(out) is [9,size(omegap)], where omegap is omega squeezed along the
% direction dim.

if nargin<2 || isempty(dim)
    dim = 1;
end

dimOmegaIn = size(omega);
oo = setdiff(1:numel(dimOmegaIn),dim);
omega = permute(omega,[dim,oo]);

ret = sqrt(sum(omega.^2,1));
PA = bsxfun(@rdivide,omega,ret);% enforce unitary rotation axis

K = cat(1,zeros(1,prod(dimOmegaIn(oo))),PA(3,:),-PA(2,:),-PA(3,:),zeros(1,prod(dimOmegaIn(oo))),PA(1,:),PA(2,:),-PA(1,:),zeros(1,prod(dimOmegaIn(oo))));
KK = cat(1,PA(1,:).^2-1,PA(1,:).*PA(2,:),PA(1,:).*PA(3,:),PA(1,:).*PA(2,:),PA(2,:).^2-1,PA(2,:).*PA(3,:),PA(1,:).*PA(3,:),PA(2,:).*PA(3,:),PA(3,:).^2-1);

inds = find(ret == 0);

out = reshape(bsxfun(@plus,bsxfun(@times,K,sin(ret(1,:))) + bsxfun(@times,KK,1-cos(ret(1,:))),[1,0,0,0,1,0,0,0,1]'),[9,dimOmegaIn(oo)]);

if ~isempty(inds)
    out(:,inds) = repmat([1;0;0;0;1;0;0;0;1],[1,numel(inds)]);% make identity matrix
end