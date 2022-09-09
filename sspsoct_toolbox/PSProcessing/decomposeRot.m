function [AOut,varargout] = decomposeRot(Min,exact)
% [AOut,varargout] = decomposeRot(Min,exact) 
% Decomposes the rotation matrix Min, in linear format, with 9 components
% along the first dimension, into the corresponding retardation vector
% Aout;
% If exact is set to true, the normalization of the optic axis is enforced.
% This is generally required, especially if the retardation values reach 0
% or pi in some locations.
% Additional management of the case when ret appraoches pi; in that case,
% the standard paloc becomse ill-defined. Instead, we extract the optic 
% axis from the diagonal of the matrix. This leaves a sign uncertainty, 
% which is resolved by leveraring the fact that this vector is an
% eigenvector and if applied with correct signs to Min is unchanged.

if nargin<2
    exact = true;
end

% This should provide the correct retardation for all possible retardations
aaa = sum(Min([1,5,9],:,:,:),1)/2-1/2;
aaa(aaa>1)=1;
aaa(aaa<-1)=-1;
retloc = real(acos(aaa));
if exact
    paloc = cat(1,(Min(6,:,:,:)-Min(8,:,:,:)),(Min(7,:,:,:)-Min(3,:,:,:)),(Min(2,:,:,:)-Min(4,:,:,:)));

    % find cases where ret approaches pi, resultin in paloc being ill-defined, because sin(ret) = 0
    mm = (pi-retloc)<eps^.125;
    if sum(mm(:))>0
        % M = I + P*sin(ret) + P*P*(1-cos(ret)); the diagonal of P*P
        % reveals the square of the pa components
        pa1 = sqrt(max((Min(1,mm)-cos(retloc(1,mm)))./(1-cos(retloc(1,mm))),0));
        pa2 = sqrt(max((Min(5,mm)-cos(retloc(1,mm)))./(1-cos(retloc(1,mm))),0));
        pa3 = sqrt(max((Min(9,mm)-cos(retloc(1,mm)))./(1-cos(retloc(1,mm))),0));

        % (M-I)*pa = pa
        temp1 = (Min(1:3,mm) - [1;0;0]).*pa1;
        temp2 = (Min(4:6,mm) - [0;1;0]).*pa2;
        temp3 = (Min(7:9,mm) - [0;0;1]).*pa3;
        
        % find which sign combination of temp is zero (or closest to zero)
        [~,mp] = min(cat(1,sum((temp1 + temp2 + temp3).^2,1),sum((-temp1 + temp2 + temp3).^2,1),sum((temp1 - temp2 + temp3).^2,1),sum((temp1 + temp2 - temp3).^2,1)),[],1);
        ss1 = 2*(mp~=2)-1;
        ss2 = 2*(mp~=3)-1;
        ss3 = 2*(mp~=4)-1;
        patemp = cat(1,pa1.*ss1,pa2.*ss2,pa3.*ss3);
        
        % there remains an overall sign ambiguity; figure out if + or - patemp matches better to paloc
        [~,mp] = min(cat(1,sum((paloc(:,mm) - patemp).^2,1),sum((paloc(:,mm) + patemp).^2,1)),[],1);
        patemp(:,mp==2) = -patemp(:,mp==2);
        
        paloc(repmat(mm,[3,1])) = patemp;
    end
    paloc = bsxfun(@rdivide,paloc,sqrt(sum(paloc.^2,1)));
    paloc(repmat(retloc==0,[3,1])) = 0;
else
    paloc = bsxfun(@times,1/2./(sin(retloc)+1e-9),cat(1,(Min(6,:,:,:)-Min(8,:,:,:)),(Min(7,:,:,:)-Min(3,:,:,:)),(Min(2,:,:,:)-Min(4,:,:,:))));
end

AOut = bsxfun(@times,paloc,retloc);

if nargout>1
    varargout{1} = retloc;
end
if nargout>2
    varargout{2} = paloc;
end
