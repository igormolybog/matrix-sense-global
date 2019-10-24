function d = LMI_P(x,z,S,W)
% implements the function LMI_P(x,z;T) from the NIPS2019 paper
% -------------------------------------------------------------------------
% x, z are matrices (n times r)
% T is a matrix
assert(all(size(x) == size(z)), 'the size of x and z dont match')
[n, r] = dims(x);
mate = x*x.' - z*z.';
e = mate(:);
if norm(e) == 0
    d = 1;
    'WARNING: x*x.T=z*z.T obtained'
    return
end
operatorX = get_X_operator(x);
smallP = P(x, z);
bigP = kron(smallP,smallP);
cvx_begin sdp quiet
    variable delta
    variable H(n*n,n*n) symmetric
    variable M(r*n,r*n) symmetric
    minimize(delta)
    operatorX.'*H*e == 0;
    (kron(eye(r),reshape(H*e,n,n)+reshape(H*e,n,n)).')+ operatorX.'*H*operatorX == M;
    M >= 0;
    S.*H == 0;
    W*H(:) == 0;
    (1 + delta)*(bigP.'*bigP) - bigP.'*H*bigP  >= 0;
    bigP.'*H*bigP - (1 - delta)*(bigP.'*bigP)  >= 0;
cvx_end
d = delta;
end

% %% Richard-like implementation
% function d = LMI_P(x,z,S,W)
%     assert(all(size(x) == size(z)), 'the size of x and z dont match')
%     [n, r] = dims(x);
%     mate = x*x.' - z*z.';
%     e2 = mate(:);
%     X2 = get_X_operator(x);
%     smallP = P(x, z);
%     bigP = kron(smallP,smallP);
%     
%     H = sdpvar(n^2,n^2);
%     SOC = sdpvar(n*r,n*r);
%     beta = sdpvar;
%     constr = [X2'*H*e2 == 0; % FOC
%             2*kron(eye(r), reshape(H*e2,n,n)) + X2'*H*X2 == SOC; SOC>=0; % SOC
%             S.*H == 0;
%             W*H(:) == 0;
%             (1 - beta)*(bigP.'*bigP) <= bigP.'*H*bigP <= (1 + beta)*(bigP.'*bigP)]; % Condition number
%     dat = optimize(constr,-beta,sdpsettings('solver','mosek','verbose',1));
%     if dat.problem == 0
%         d = beta;
%     else
%         d = inf;
%     end
% end