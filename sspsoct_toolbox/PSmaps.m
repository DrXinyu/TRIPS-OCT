function [ret1,oauim,dop] = PSmaps(S1,S2,fwx,dz,CalStru)

%% filtering
I1 = sqrt(dot(S1,S1,4));
I2 = sqrt(dot(S2,S2,4));

nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
h = exp(-nx.^2);
h = h/sum(h(:));

S1f = imfilter(S1,h,'circular');
S2f = imfilter(S2,h,'circular');
I1f = imfilter(I1,h,'circular');
I2f = imfilter(I2,h,'circular');


% final normalization
QUVf1 = sqrt(dot(S1f,S1f,4));
QUVf2 = sqrt(dot(S2f,S2f,4));

dop1 = QUVf1./I1f;
dop2 = QUVf2./I2f;


S1f = S1f./repmat(QUVf1,[1,1,1,3]);
S2f = S2f./repmat(QUVf2,[1,1,1,3]);


% dop = ((dop1)+(dop2))<1.1;
% S1f(repmat(dop,[1,1,1,3])) = 0;
% S2f(repmat(dop,[1,1,1,3])) = 0;
% 
% force the two Stokes vectors to be orthogonal, which is equivalent to the
% lsq solution

%construct orthonormal tripod for these data points
na = S1f;
nb = cross(S1f,S2f,4);

S1f = na./repmat(max(sqrt(dot(na,na,4)),1e-9),[1,1,1,3]);
S2f = nb./repmat(max(sqrt(dot(nb,nb,4)),1e-9),[1,1,1,3]);




% local birefringence analysis
S1plus = circshift(S1f,-dz,1);
S1minus = circshift(S1f,dz,1); 
S2plus = circshift(S2f,-dz,1);
S2minus = circshift(S2f,dz,1);   

% simple cross product of difference vectors to find PA
PA = cross(S1minus-S1plus,S2minus-S2plus,4);
PA = PA./repmat(max(sqrt(dot(PA,PA,4)),1e-9),[1,1,1,3]);

temp = dot(S1plus,PA,4).^2;% subtract tiny number to avoid division by zero
%in the next line
%retSinW = real(acos((dot(S1plus,S1minus,4)-temp)./(1-temp + 1e-9)));
retSinW = real(acos(complex(dot(S1plus,S1minus,4)-temp)./(1-temp + 1e-9)));

pm = sign((1-dot(S1plus,S1minus,4)).*(dot(S1plus-S1minus,S2plus+S2minus,4)));
PA = PA.*repmat(pm,[1,1,1,3]);

taucorr = gather(PA.*retSinW);

%% weight for spectral alignment
dim = size(taucorr);
mid = ceil(dim(3)/2);
weight = ones(size(retSinW)).*dop1.*dop2;%.*(QUVf1+QUVf2);%log((QUVf2+QUVf1)./sum(QUVf1+QUVf2,3));
weight = gather(weight);
avgWeg=weight;
B = squeeze(taucorr(:,:,mid,:));
for wind = [(1:mid-1),(mid+1:dim(3))]
    ww = weight(:,:,wind).*weight(:,:,mid);
    A = squeeze(bsxfun(@times,taucorr(:,:,wind,:),ww));
    
    C(1,1) = sum(sum(A(:,:,1).*B(:,:,1)));
    C(2,1) = sum(sum(A(:,:,1).*B(:,:,2)));
    C(3,1) = sum(sum(A(:,:,1).*B(:,:,3)));
    
    C(1,2) = sum(sum(A(:,:,2).*B(:,:,1)));
    C(2,2) = sum(sum(A(:,:,2).*B(:,:,2)));
    C(3,2) = sum(sum(A(:,:,2).*B(:,:,3)));
    
    C(1,3) = sum(sum(A(:,:,3).*B(:,:,1)));
    C(2,3) = sum(sum(A(:,:,3).*B(:,:,2)));
    C(3,3) = sum(sum(A(:,:,3).*B(:,:,3)));
    
    R = reshape(euclideanRotation(C(:)),[3,3]);
    
    temp1 = R(1,1)*taucorr(:,:,wind,1) + R(1,2)*taucorr(:,:,wind,2) + R(1,3)*taucorr(:,:,wind,3);
    temp2 = R(2,1)*taucorr(:,:,wind,1) + R(2,2)*taucorr(:,:,wind,2) + R(2,3)*taucorr(:,:,wind,3);
    temp3 = R(3,1)*taucorr(:,:,wind,1) + R(3,2)*taucorr(:,:,wind,2) + R(3,3)*taucorr(:,:,wind,3);
    
    taucorr(:,:,wind,1) = temp1;
    taucorr(:,:,wind,2) = temp2;
    taucorr(:,:,wind,3) = temp3;
    rcorrglobal(:,wind) = decomposeRot(R(:));
end
avgWeg(avgWeg<=0) = 1e-9;
taumean = squeeze(sum(bsxfun(@times,taucorr,avgWeg),3));
taumean = bsxfun(@rdivide, taumean,sum(avgWeg,3));

taumean(isnan(taumean)) = 0;
% tauErr = sqrt(sum((bsxfun(@minus,taucorr,permute(taumean,[1 2 4 3]))).^2,4));
% tauErr(repmat(mask,[1 1 3])) = 0;
% out.binError1 = sum(squeeze(sum(tauErr)));
% 
ret1 = sqrt(sum(taumean.^2,3));%*0.0239/2/dz;
unit_depth = (CalStru.whole_depth/CalStru.NFFT)*dz;
ret1 = ((ret1/(2*pi))*360)/unit_depth/2;

oaim = taumean./repmat(max(sqrt(dot(taumean,taumean,3)),1e-9),[1,1,3]);
oauim = atan2(oaim(:,:,1),oaim(:,:,2));
dop = gather(mean(dop1+dop2,3));






end

