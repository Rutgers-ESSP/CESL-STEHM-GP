function [alfa,invtraincv,s,r]=svdinv(traincv,y0)

% [alfa,invtraincv,s,r]=svdinv(traincv,y0)
%
% Do a SVD of traincv, calculate the inverse of traincv,
% then calculate alfa = invtraincv * y0.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 28 16:31:00 EST 2015
    
    [m,n] = size(traincv);
    [U,S,V] = svd(traincv,0);
    s = diag(S);
    tol = max(m,n) * eps(max(s));
    r = sum(s > tol);
    invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
    alfa = invtraincv * y0;