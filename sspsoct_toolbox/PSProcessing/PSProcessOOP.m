function out = PSProcessOOP(S1,S2,procStruct)
%PSPROCESS - computes local birefringence of the Stokes vectors S1 and S2,
%but accepts four dimensional input Stokes vectors, i.e. size(S1,4) = 4, to
%enable previous out of plane averaging of the Stokes vectors.
% 
% The third dimension is the pc multiplexed dimension. 
% Using the dop is default for defining the correction mask
% Also, using pcEst is default

% the only two mandatory arguments
fwx = procStruct.fwx;
dz = procStruct.dz;


if ~isfield(procStruct,'fwz')
    fwz = 1;
else
    fwz = procStruct.fwz;
end

% implementation of estimation of the relative rotation between the PAs of
% the subwindows
if isfield(procStruct,'pcEst')
    pcEst = procStruct.pcEst;
else
    pcEst = true; % pcEst is default
end

if isfield(procStruct,'dopTh')
    dopTh = procStruct.dopTh;
else
    dopTh = [0.6,1];
end

% gaussian lateral filtering is default
if ~isfield(procStruct,'filter')
    procStruct.filter = 'gaussian';
end


% padd Stokes vectors on both sides
prepad = floor(ceil(1.5*fwx)/2)+floor(fwx/2);
postpad = ceil(ceil(1.5*fwx)/2)+ceil(fwx/2);

S1 = cat(2,S1(:,end-prepad+1:end,:,:),S1,S1(:,1:postpad,:,:));
S2 = cat(2,S2(:,end-prepad+1:end,:,:),S2,S2(:,1:postpad,:,:));

% 
% I1 = sqrt(dot(S1,S1,4));
% I2 = sqrt(dot(S2,S2,4));
% 
% if dopuCheck
%     % pre-normalization for computation of DOPU
%     Sn1 = S1./repmat(I1,[1,1,1,3]);
%     Sn2 = S2./repmat(I2,[1,1,1,3]);
% end


%II = i1 + i2;

% filtering of S1 and S2
if isfield(procStruct,'filter') && strcmp(procStruct.filter,'mean') 
    S1 = imfilter(S1,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
    S2 = imfilter(S2,ones(fwz,fwx,fwy)/fwx/fwy/fwz,'circular');
elseif isfield(procStruct,'filter') && strcmp(procStruct.filter,'gaussian') 
    nx = (round(fwx*1.5)-1)/2;
    nz = (round(fwz*1.5)-1)/2;
    nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
    nz = linspace(-nz,nz,round(fwz*1.5))'*2*sqrt(log(2))/fwz;
    if fwz == 1
        nz = 0;
    end
    h = exp(-nz.^2)*exp(-nx.^2);
    h = h/sum(h(:));
    h = repmat(h,[1,1,fwy]);% out of plane it remains a simple mean filter
    S1 = imfilter(S1,h,'circular');
    S2 = imfilter(S2,h,'circular');
end


% final normalization
I1 = S1(:,:,:,1);
I2 = S2(:,:,:,1);
S1 = S1(:,:,:,2:4);
S2 = S2(:,:,:,2:4);
L1 = sum(S1.^2,4);
L2 = sum(S2.^2,4);

dop = sqrt(sum(L1 + L2,3)./sum(I1.^2 + I2.^2,3));

If = sqrt(I1.^2 + I2.^2);

S1 = S1./repmat(If1,[1,1,1,3]);
S2 = S2./repmat(If2,[1,1,1,3]);

% force the two Stokes vectors to be orthogonal, which is equivalent to the
% lsq solution

% construct orthonormal tripod for these data points
na = S1 + S2;
nb = S1 - S2;
S1 = na./repmat(sqrt(dot(na,na,4)),[1,1,1,3]);
S2 = nb./repmat(sqrt(dot(nb,nb,4)),[1,1,1,3]);

% local birefringence analysis
S1plus = circshift(S1,-dz);
S1minus = circshift(S1,dz); 
S2plus = circshift(S2,-dz);
S2minus = circshift(S2,dz);   

% simple cross product of difference vectors to find PA
PA = cross(S1minus-S1plus,S2minus-S2plus,4);
PA = PA./repmat(max(sqrt(dot(PA,PA,4)),1e-9),[1,1,1,3]);

temp = dot(S1plus,PA,4).^2;
retSinW = real(acos((dot(S1plus,S1minus,4)-temp)./(1-temp)))/2/dz;
pm = sign((1-dot(S1plus,S1minus,4)).*(dot(S1plus-S1minus,S2plus+S2minus,4)));
PA = PA.*repmat(pm,[1,1,1,3]);

rpamean = mean(PA.*repmat(retSinW,[1,1,1,3]),3);
rmean = sqrt(dot(rpamean,rpamean,4));
out.stdpa = 1-rmean(:,prepad+1:end-postpad)./mean(abs(retSinW(:,prepad+1:end-postpad,:)),3);    
out.rmean = rmean(:,prepad+1:end-postpad)*100/4.8*180/pi;


if pcEst
    % implementation of estimation of relative rotations between the
    % different subwindows

    Wcorr = PA.*repmat(retSinW,[1,1,1,3]);
    
    dim = size(PA);

    mask = (dop>dopTh(1)).*(dop<=dopTh(2));
     
%    mask(1:100,:) = 0;
%    mask(801:end,:) = 0;
    
    mid = ceil(size(retSinW,3)/2);
    PA = PA.*repmat(mask,[1,1,dim(3),3]);

    ref = PA(:,:,mid,:);
    C = zeros([3,3,dim(2)]);

    h = ones(1,fwx)/3;
    for wind = [(1:mid-1),(mid+1:dim(3))]
        C(1,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,1),1),h,'same');
        C(2,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,2),1),h,'same');
        C(3,1,:) = conv(sum(PA(:,:,wind,1).*ref(:,:,3),1),h,'same');

        C(1,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,1),1),h,'same');
        C(2,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,2),1),h,'same');
        C(3,2,:) = conv(sum(PA(:,:,wind,2).*ref(:,:,3),1),h,'same');

        C(1,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,1),1),h,'same');
        C(2,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,2),1),h,'same');
        C(3,3,:) = conv(sum(PA(:,:,wind,3).*ref(:,:,3),1),h,'same');

        R = reshape(euclideanRotation(reshape(C,[9,size(C,3)])),[3,3,size(C,3)]);

        temp1 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(1,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(1,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(1,3,:),1),[dim(1),1]);
        temp2 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(2,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(2,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(2,3,:),1),[dim(1),1]);
        temp3 = Wcorr(:,:,wind,1).*repmat(shiftdim(R(3,1,:),1),[dim(1),1]) + Wcorr(:,:,wind,2).*repmat(shiftdim(R(3,2,:),1),[dim(1),1]) + Wcorr(:,:,wind,3).*repmat(shiftdim(R(3,3,:),1),[dim(1),1]);

        Wcorr(:,:,wind,1) = temp1;
        Wcorr(:,:,wind,2) = temp2;
        Wcorr(:,:,wind,3) = temp3;
    end
    
    rpamean = mean(Wcorr,3);
    rmeancorr = sqrt(dot(rpamean,rpamean,4));

end

if pcEst
    out.rmeancorr = rmeancorr(:,prepad+1:end-postpad)*100/4.8*180/pi;
    out.PAcorr = rpamean(:,prepad+1:end-postpad,:,:)./repmat(rmeancorr(:,prepad+1:end-postpad),[1,1,1,3]);
    %    out.errCorr = errCorr;
end

%out.I = squeeze(II(:,prepad+1:end-postpad,:));
out.If = squeeze(mean(I1(:,prepad+1:end-postpad,:)+I1(:,prepad+1:end-postpad,:),3));
out.ret = squeeze(retSinW(:,prepad+1:end-postpad,:))*100/4.8*180/pi;
out.PA = squeeze(PA(:,prepad+1:end-postpad,:,:));
    
out.dop = dop(:,prepad+1:end-postpad,:);



