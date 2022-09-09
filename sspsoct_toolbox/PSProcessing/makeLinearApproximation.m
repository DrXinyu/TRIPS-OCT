function Mout = makeLinearApproximation(Min)
% Mout = makeLinearApproximation(Min) transforms the deterministic input 
% Mueller matrix Min into a purely linear approximation in the least-square
% sense, which corresponds to the H-filtering of the Min + G*Min.'*G*.

dim = size(Min);
if dim(1) == 4 && dim(2) == 4
    if numel(dim)>2
        Min = reshape(Min,cat(2,16,dim(3:end)));
    else
        dim = [16,1];
        Min = Min(:);
    end
end

% convert M to H matrix; Since the input is assumed to be deterministic, H
% has rank 1 and the first column is sufficient to deduce the unerlying
% Jones matrix.
A = 1/sqrt(2)*[1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0];
At = A';

H = MatrixMultiply(At(:),MatrixMultiply(Min,A(:)));
H([2,4,6,8,9,11,13,15],:) = H([9,11,13,15,2,4,6,8],:);
J = H([1;3;2;4],:)./sqrt(H(1,:));
% J is the linearized Jones matrix, however in row order, i.e. J =
% [j11;j12;j21;j22].

% let's construct an orthogonal basis from J and its transpose
w1 = J./sqrt(sum(abs(J).^2,1));
w2 = w1([1;3;2;4],:) - sum(conj(w1).*w1([1;3;2;4],:),1).*w1;
w2 = w2./sqrt(sum(abs(w2).^2,1));

% manage the case where w1 is already linear by setting w2 to zero
w2(:,abs(w1(2,:) - w1(3,:))<eps) = 0;

% We are now looking for a linear combination of w1 and w2 that defines the
% direction that corresponds that achieves the maximum projection of J and
% its transpose. This corresponds to the first principle component, or the
% eigenvector corresponding to the highest eigenvalue.
% wopt = [w1 w2]*[a;b], where a.^2 + b.^2 =1, a is real and b is complex-
% valued (an overall phase term is irrelevant). We thus look for 
% argmin abs([J,Jtranspose]'*[w1 w2]*[a;b]) =
% [a;b]'*[w1,w2]'*[J,Jt]*[J,Jt]'*[w1,w2]*[a;b] = [a;b]'*W*[a;b]
% Parameterizing a = cos(alpha) and b = sin(alpha)*exp(1i*beta), we find
% analytical solutions for alpha and beta:

W = cat(1,sqrt(sum(abs(J).^2)),sum(conj(J([1;3;2;4],:)).*w1,1),sum(conj(J).*w2),sum(conj(J([1;3;2;4],:)).*w2));
W = MatrixMultiply(conj(W([1;3;2;4],:)),W);

beta = atan2(imag(W(2,:)),real(W(2,:)));
alpha = atan2(2*(real(W(2,:)).*cos(beta) + imag(W(2,:)).*sin(beta)),real(W(1,:)-W(4,:)))/2;

Jout = cos(alpha).*w1 +  sin(alpha).*exp(1i*beta).*w2;

% figure out corresponding eigenvalue to preserve scaling of the Mueller
% matrix
ev = sqrt(sum(conj(Jout).*(H(1:4,:).*Jout(1,:) + H(5:8,:).*Jout(2,:) + H(9:12,:).*Jout(3,:) + H(13:16,:).*Jout(4,:)),1));
Jout = Jout.*abs(ev);

% Convert to Mueller
Mout = Jones2Mueller(Jout([1;3;2;4],:));
Mout = reshape(Mout,dim);

