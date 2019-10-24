% Initial setup
addpath('util/');

T = 72000;
n = 8;
r = 3;


%%
% Construction of H for the ellipsoid norm structure

blocks = cell(n);
for i = 1:n
    blocks{i} = i*ones(n);
end
H = blkdiag(blocks{:});
[S,W] = get_structure(H);


%%
% Tossing

'points that return bound < 1:'

for toss = 1:T
    x = randi([-1 1], n, r);
    y = randi([-1 1], n, r);
    delta = LMI(x, z, S, W);
    if delta < 1
        toss, delta, x, y
        break
    end
end

'.. and that is it'