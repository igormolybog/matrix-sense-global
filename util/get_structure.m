function [S,W] = get_structure(H_ish)
% returns the matrix representation of the components of the structure 
% operator introduced in Section 3
S = ones(size(H_ish)) - spones(H_ish);
W = get_identity_structure(H_ish);
end