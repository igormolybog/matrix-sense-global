% Initial setup

n = 8;
r = 2;
T = 1000;

%%
% tossing

for toss = 1:T
    x = randi([-1 1], n, r);
    y = randi([-1 1], n, r);
    a = counterixample_if_greater_than_zero(x, y);
    if a > 0.1
        toss, a, x, y
        break
    end
end

%%
% Solving the LMI without a constraint on the hessian

function a = counterixample_if_greater_than_zero(x, z)
    n = size(x,1);
    r = size(x,2);
    D = z*z'-x*x';
    cvx_begin sdp quiet
        variable a
        variable M(n,n) symmetric
        minimize(-a)
        M-a*eye(n) >= 0;
        (D*M + M*D)*z == 0;
        trace(D*M*D) >= 0.1;
    cvx_end
    cvx_status
end