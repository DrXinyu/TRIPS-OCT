function ball = decomposeBallLens(surf)
% ball = decomposeBallLens(surf) takes the optic axis evolution of
% the signal from the ball lens (or the inner sheath interface) to 
% decompose it into M.'*M where M = L2*V*L1 (Jones), L1 and L2 are linear
% retarders, and V is a circular retarder.
% vcomp is a boolean, if true, all angle offsets occurring in V are 
% compensated, if false, only the linear part is.
% Under development; uses trace and rigorous LSQ formulation


% input:
% surf : the correctly unwrapped retardation vector of the ball lens signal
MSurf = makeRot(surf);


% optimize L1 to reduce the overall trace error, keep Lq2 zero and Lq3
% constant
fun = @(x) findLsTrace(x,MSurf);
start = double(mean(surf,2)/2);
xopt = fminsearch(fun,start);
[err,OA] = fun(xopt);

% check if the resulting angle of the isolated ball lens signal increases 
% add pi to xopt otherwise to result in a flip (this is independent of
% surf, as it uses directly the ambiguity-free matrix Msurf) for
% computation).
phi = unwrap(atan2(OA(2,:),OA(1,:)));
ball.L1flip = false;
if mean(diff(phi))<0% negative rotation
    xopt(1:2) = xopt(1:2) - pi*bsxfun(@rdivide,xopt(1:2),norm(xopt(1:2)));
    [err,OA] = fun(xopt);
    ball.L1flip = true;
end
OARet = squeeze(sqrt(sum(OA.^2,1)));

[~,~,oamodel,Vcorr,mret,wrap] = fun(xopt);
% check of oamodel aligned with surf signal, or if the alternate square
% root solution should be computed.
ball.rootflip = false;
if wrap
    surf = surf - 2*pi*bsxfun(@rdivide,surf,sqrt(sum(surf.^2,1)));
    ball.rootflip = true;
end
MVcorr = makeRot(cat(1,zeros(2,size(surf,2)),Vcorr));
Ballroot = makeRot(-surf/2);
Mcorr = MatrixMultiply(Ballroot,MVcorr);

ball.Mcorr = Mcorr;
ball.OA = OA; % unwrapped central elements
ball.RL1 = xopt(1:2); % Linear elements L1
ball.oaret = OARet;
ball.Vcomp = Vcorr; 
ball.surf = surf;
ball.oamodel = oamodel;
ball.err = err;
ball.mret = mret;


function [err,oacenter,oamodel,Vcorr,mret,wrapFlag] = findLsTrace(params,Mball)
% [err,oacenter,oamodel,Vcorr,mret,wrapFlag] = findLsTrace(params,Mball) 
% computes the error of the lsq trace minimization problem. oacenter is the
% optic axis of the ball lens signal after correcting with a linear element
% L1. oamodel is the optic axis of the modeled double-pass transmission, 
% Vcorr, the rotation required to apply to the signal after compensating it
% with the inverse square root of the ball lens signal.
% The model corresponds to L1*W*L2*L2*Winv*L1

N = size(Mball,2);

% This first part simply compensates Mball with L1, and computes the best
% L2, which always aligns to Mball by the actcion of W.
% construct candidate L1 matrix
L1inv = makeRot(-[params(1);params(2);0]);
Mcomp = MatrixMultiply(L1inv,MatrixMultiply(Mball,L1inv));

% compute unwrapped ret and then compute the trace minimum
oacenter = decomposeRot(Mcomp); % unwrapping is necessary, as this is an artifact of the varying optic axis orientation; 
ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(oacenter,[],2).^2,1))>pi,2),2)),[3,1])>0;
delta = 2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
oacenter(ss) = oacenter(ss) - delta(ss);
ret = sqrt(sum(oacenter.^2,1));
if mean(ret)>pi % the modulo 2*pi solutions do not change the sense of orientation, and simply switch the sign of mret
    oacenter = oacenter-2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
    ret = sqrt(sum(oacenter.^2,1));
end
mret = atan2(mean(sin(real(ret))),mean(cos(real(ret))));
err = 1-mean(cos(ret-mret));


if nargout>2
    % This assumes that the V direction was adjusted to be positive

    % Compute angle of rotation matrix W
    angleW = atan2(oacenter(2,:),oacenter(1,:));
    Winv = makeRot(cat(1,zeros(2,N),-angleW));
    L2 = makeRot([mret/2;zeros(2,numel(mret))]);
    L1 = makeRot([params(1);params(2);0]);
    Mfw = MatrixMultiply(L2,MatrixMultiply(Winv,L1));
    oamodel = decomposeRot(MatrixMultiply(bsxfun(@times,Mfw([1,4,7,2,5,8,3,6,9]',:),[1;1;-1;1;1;-1;-1;-1;1]),Mfw));
    % unwrapping to avoid inconsistencies
    ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(oamodel,[],2).^2,1))>pi,2),2)),[3,1])>0;
    delta = 2*pi*bsxfun(@rdivide,oamodel,sqrt(sum(oamodel.^2,1)));
    oamodel(ss) = oamodel(ss) - delta(ss);
    if mean(sqrt(sum(oamodel.^2,1)))>pi % this is only necessary for consistency with surf, which was unwrapped in the same way
        oamodel = oamodel-2*pi*bsxfun(@rdivide,oamodel,sqrt(sum(oamodel.^2,1)));
    end
    % the square root allows to solutions; one of them will result in a
    % flipping of the recovered optic axes. We have to figure out which one
    % matches the computed model, and then use the matching root to
    % construct the correction matrix.
    Msqrinv = makeRot(-oamodel/2);
    V1 = decomposeRot(MatrixMultiply(Mfw,Msqrinv));
    oamodel2 = oamodel-2*pi*bsxfun(@rdivide,oamodel,sqrt(sum(oamodel.^2,1)));
    Msqrinv = makeRot(-oamodel2/2);
    V2 = decomposeRot(MatrixMultiply(Mfw,Msqrinv));
    % the correct solution has only components along V
    if mean(abs(V1(3,:)./sqrt(sum(V1.^2,1))),'omitnan') > mean(abs(V2(3,:)./sqrt(sum(V2.^2,1))),'omitnan')
        Vcorr = -unwrap(V1(3,:));
        wrapFlag = false;
    else
        Vcorr = -unwrap(V2(3,:));
        wrapFlag = true;
    end
    
end

