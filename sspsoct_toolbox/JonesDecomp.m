function [r,d] = JonesDecomp(J,polar,pureDiatt)
% JonesDecomp computes the retardation and diattenuation from a Jones
% matrix, provided as input argument ( the matrix logarithm is part of the
% function, using the 'concurrent' decomposition.
% If called with the additioanl argument 'polar', it performs the polar
% decomposition, decomposing the Jones matrix into a sequence of a
% diattenuation and a retardation matrix, J = Jr*Jd.
% Different approach without logm
% J is either 2x2
% or 4 x whatever, where the 4 elements are [J11;J21;J12;J22]

dim = size(J);

if nargin<2
    polar = false;
end
if nargin<3
    pureDiatt = false;
end

if ~polar
    if dim(1) == 2 && dim(2) == 2
        J = J/sqrt(det(J));
        q = cat(1,(J(1,1)-J(2,2)),(J(2,1)+J(1,2)),(-1i*J(2,1)+1i*J(1,2)))/2;
        tr = trace(J)/2;
        c = acosh(tr);
        csin = c/sinh(c);
        csin(c==0) = 1;
        f = 2*q*csin;
        r = -imag(f);
        d = real(f);

    elseif dim(1)==4
        detJ = sqrt(J(1,:).*J(4,:)-J(2,:).*J(3,:));
        J = bsxfun(@rdivide,J(:,:),detJ);
        q = cat(1,(J(1,:)-J(4,:)),(J(2,:)+J(3,:)),(-1i*J(2,:)+1i*J(3,:)))/2;
        tr = (J(1,:) + J(4,:))/2;
        c = acosh(tr);
        csin = c./sinh(c);
        csin(c==0) = 1;
        f = 2*bsxfun(@times,q,csin);
%        f = 2*bsxfun(@times,q,c./sinh(c));
        r = reshape(-imag(f),[3,dim(2:end)]);
        d = reshape(real(f),[3,dim(2:end)]);
    end
else% polar decomposition
    if dim(1) == 2 && dim(2) == 2
        J = J/sqrt(det(J));
        Jd = J'*J;
        [~,d] = JonesDecomp(Jd);
        d = d/2;
        r = JonesDecomp(J*makeJones(zeros(size(d)),-d));
    elseif dim(1)==4
        detJ = sqrt(J(1,:).*J(4,:)-J(2,:).*J(3,:));
        J = bsxfun(@rdivide,J(:,:),detJ);
        [~,d] = JonesDecomp(MatrixMultiply(conj(J([1;3;2;4],:)),J));
        d = d/2;
        r = JonesDecomp(MatrixMultiply(J,makeJones(zeros(size(d)),-d)));
        r = reshape(r,[3,dim(2:end)]);
        d = reshape(d,[3,dim(2:end)]);
    end
end
