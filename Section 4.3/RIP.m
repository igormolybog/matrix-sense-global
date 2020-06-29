function [RIP, xmi, zmi, xma, zma] = RIP(H, n_samples)
% Estimates the RIP constant that corresponds to the kernel operator H
ma = 0;
mi = inf;

nb = sqrt(size(H,2));
for i = 1:n_samples
   x = normrnd(1,0.1,nb,1);
   z = -normrnd(1,0.1,nb,1);
   sample = sampleRI(H, [x,z]);
   if sample > ma
    ma = max(ma, sample);
    xma = x;
    zma = z;
   end
   if sample < mi
    mi = min(mi, sample);
    xmi = x;
    zmi = z;
   end
end

RIP = (ma/(mi+ma))*2-1; %RIP constant of A

end

function sample = sampleRI(H, X)
    XX = X*X.';
    sample = (XX(:)'*H*XX(:)/norm(X,'F'))^2;
end
