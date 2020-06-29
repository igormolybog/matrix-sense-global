function H = buildH(y)
% y --- vector (COLUMN) y_1, .. y_{n-1}

if size(y,1)<size(y,2)
    y = y';
end
assert(size(y,2) == 1);
n = size(y,1)+1;
H = sparse(n*n, n*n);
for i = 1:n-1
H = H + y(i)^2*vec(E(i,i,n)-E(n,i,n))*vec(E(i,i,n)-E(n,i,n))';
end

y_n = sum(y);
a_n = sparse(n,n);
for j = 1:n-1
    a_n = a_n - y(j)*E(j,n,n); 
end
a_n = a_n + y_n*E(n,n,n);

H = H + vec(a_n)*vec(a_n)';
end

function E = E(i,j,n)
E = sparse(n,n);
E(i,j) = 1;
end

function vec = vec(Matrix)
vec = Matrix(:);
end

% function unvec = unvec(Vector)
% unvec = reshape
% end