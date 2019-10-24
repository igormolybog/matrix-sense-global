function W = get_identity_structure(H_ish)
% returns the matrix representation of a linear operator (acting on vectorized version of H, although 
% the function receives a non-vectorized symmetric matrix H_ish) having null-space with the property:
% if two components of H_ish has the same values, then the corresponding two
% components of any vector in the null-space must have same value as well
% -----------------------------------------------------------------------
% Example: 
% [ 1   0   0.5]
% [ 0  0.5   1 ] = H_ish
% [0.5  1    -1]
%  
% vectorized verstion of
% [ 3   3   100]
% [ 3  100   1 ]   is in the null-space of W
% [100  1    0 ]
% 
% vectorized version of
% [ 3   3   100]
% [ 3   1    1 ]   is not in the null-space of W
% [100  1    0 ]


    triuH_ish = triu(H_ish);
    value = unique(triuH_ish);
    veclen = size(H_ish,1)*size(H_ish,1);
    W = sparse(0, veclen);
    for i = 1:size(value, 1)
        if all(value(i) == 0)
            continue
        end
        comps = find(triuH_ish==value(i));
        for j = 2:size(comps, 1)
            newline = sparse([1 1], [comps(1), comps(j)], [1 -1], 1, veclen);
            W = [W; newline];
        end
    end
    if size(W,1) == 0 
        W = 0;
    end
end