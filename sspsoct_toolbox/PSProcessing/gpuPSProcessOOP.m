function [retout,dopout,varargout] = gpuPSProcessOOP(S1,S2,fwx,dz,Nw,prepad,postpad,resz)
% Reconstructs the local retardation and DOP measures with the spectral
% binning algorithm on the GPU. The input has to be already padded for
% circular or other lateral filtering; i.e. input arrays are prepad+postpad
% wider than the output
% To enable out of plane averaging, the Stokes vectors are previously
% averaged, and passed in with all four components, i.e. size(S1,4) = 4
% The optional third output argument returns the filtered intensity (summed
% over all spectral bins).

% All arrays are 2D. The ordering
% of the input Stokes vectors is Nk,NAlines,Nw,[S,Q,U,V]
%
% Arguments:
% S1,S2; Stokes vectors
% fwx: lateral filter fwhm
% dz: half width of axial offset for derivation of the local retardation
% Nw: number of windows present
% prepad, postpad, number of padding to be removed of output
% resz: scaling of axial pixel, in um/pixel (default: 4.8)

if nargin<8
    resz = 4.8;
end

% threshold on dop applied for PA-correction
dopTh = [0.6,1];

% construction of the lateral filter
nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
h = exp(-nx.^2);
h = h(:)/sum(h(:));

% let's permute the dimensions, as we are working along rows here
S1 = permute(S1,[2,1,3:numel(size(S1))]);
S2 = permute(S2,[2,1,3:numel(size(S2))]);

dimIn = size(S1);

% filtering with h
S1 = conv2(S1(:,:),h,'same');
S2 = conv2(S2(:,:),h,'same');

S1 = reshape(S1,dimIn);
S2 = reshape(S2,dimIn);

% Euclidian length of Q,U,V
L1 = (sum(S1(:,:,:,2:4).^2,4));
L2 = (sum(S2(:,:,:,2:4).^2,4));

If = mean(S1(:,:,:,1).^2 + S2(:,:,:,1).^2,3); % keep filtered intensity as optional output argument

% computation of DOP
dop = sqrt(mean(L1 + L2,3)./If);
If = sqrt(If);

% final normalization
S1 = bsxfun(@rdivide,S1(:,:,:,2:4),sqrt(L1));
S2 = bsxfun(@rdivide,S2(:,:,:,2:4),sqrt(L2));

clear L1 L2;

% force the two Stokes vectors to be orthogonal, which is equivalent to the
% lsq solution

% construct orthonormal tripod for these data points; efficient use of the
% arrays S1, S2 to avoid intialization of additional variables
S1 = S1 + S2;
S2 = S1 - 2*S2;

nna = sqrt(sum(S1.^2,4));
nnb = sqrt(sum(S2.^2,4));

S1 = bsxfun(@rdivide,S1,nna);
S2 = bsxfun(@rdivide,S2,nnb);

%S1(isnan(S1)) = 0;
%S2(isnan(S2)) = 0;

% local birefringence analysis
S1plus = circshift(S1,[0,-dz]);
S1 = circshift(S1,[0,dz]); 
S2plus = circshift(S2,[0,-dz]);
S2 = circshift(S2,[0,dz]);   

%PA = cross(S1-S1plus,S2-S2plus,4);
PA = gpuArray(zeros(size(S1),'like',S1));
PA(:,:,:,1) = bsxfun(@times,S1(:,:,:,2)-S1plus(:,:,:,2),S2(:,:,:,3)-S2plus(:,:,:,3)) - bsxfun(@times,S1(:,:,:,3)-S1plus(:,:,:,3),S2(:,:,:,2)-S2plus(:,:,:,2));
PA(:,:,:,2) = bsxfun(@times,S1(:,:,:,3)-S1plus(:,:,:,3),S2(:,:,:,1)-S2plus(:,:,:,1)) - bsxfun(@times,S1(:,:,:,1)-S1plus(:,:,:,1),S2(:,:,:,3)-S2plus(:,:,:,3));
PA(:,:,:,3) = bsxfun(@times,S1(:,:,:,1)-S1plus(:,:,:,1),S2(:,:,:,2)-S2plus(:,:,:,2)) - bsxfun(@times,S1(:,:,:,2)-S1plus(:,:,:,2),S2(:,:,:,1)-S2plus(:,:,:,1));

PA = bsxfun(@rdivide,PA,max(sqrt(sum(PA.^2,4)),1e-15));

temp = dot(S1plus,PA,4).^2;
temp = min(max((dot(S1plus,S1,4)-temp)./(1-temp),-1),1);
retSinW = acos(temp)/2/dz;
pm = sign((1-dot(S1plus,S1,4)).*(dot(S1plus-S1,S2plus+S2,4)));
PA = bsxfun(@times,PA,pm);

PAW = bsxfun(@times,PA,retSinW);

clear S1 S1plus S2 S2plus;
% mask of dop 
dim = size(PA);
mask = (dop>dopTh(1)).*(dop<=dopTh(2));
mid = ceil(Nw/2);

PA = bsxfun(@times,PA,mask);

% optic axis of central bin
ref = PA(:,:,mid,:);

h = ones(fwx,1)/3;
for wind = [(1:mid-1),(mid+1:Nw)]
    C = gpuArray(zeros([3,3,dim(1)],'single'));

     C(1,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,1),2),h,'same');
     C(2,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,2),2),h,'same');
     C(3,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,3),2),h,'same');
     
     C(1,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,1),2),h,'same');
     C(2,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,2),2),h,'same');
     C(3,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,3),2),h,'same');
 
     C(1,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,1),2),h,'same');
     C(2,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,2),2),h,'same');
     C(3,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,3),2),h,'same');
   
    C = gather(C);
    R = reshape(euclideanRotation(reshape(C,[9,size(C,3)])),[3,3,size(C,3)]);

    R = permute(R,[3,1,2]);
    temp1 = bsxfun(@times,PAW(:,:,wind,1),R(:,1,1)) + bsxfun(@times,PAW(:,:,wind,2),R(:,1,2)) + bsxfun(@times,PAW(:,:,wind,3),R(:,1,3));
    temp2 = bsxfun(@times,PAW(:,:,wind,1),R(:,2,1)) + bsxfun(@times,PAW(:,:,wind,2),R(:,2,2)) + bsxfun(@times,PAW(:,:,wind,3),R(:,2,3));
    temp3 = bsxfun(@times,PAW(:,:,wind,1),R(:,3,1)) + bsxfun(@times,PAW(:,:,wind,2),R(:,3,2)) + bsxfun(@times,PAW(:,:,wind,3),R(:,3,3));

    PAW(:,:,wind,1) = temp1;
    PAW(:,:,wind,2) = temp2;
    PAW(:,:,wind,3) = temp3;
end
    
%rpamean = bsxfun(@rdivide,sum(PAW,3),Nw);
rpamean = sum(PAW,3)/Nw;
rmean = sqrt(dot(rpamean,rpamean,4));

%rpamean = permute(rpamean,[2,1,4,3]);
rmean = permute(rmean,[2,1]);
dop = permute(dop,[2,1]);
If = permute(If,[2,1]);

retout = gather(rmean(:,prepad+1:end-postpad)*100/resz*180/pi);
dopout = gather(dop(:,prepad+1:end-postpad));
if nargout>2
    varargout{1} = gather(If(:,prepad+1:end-postpad));
end