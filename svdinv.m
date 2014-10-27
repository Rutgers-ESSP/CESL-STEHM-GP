
function [alfa,invtraincv,s,r]=svdinv(traincv,y0)
	[m,n] = size(traincv);
	[U,S,V] = svd(traincv,0);
	s = diag(S);
	tol = max(m,n) * eps(max(s));
	r = sum(s > tol);
	invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
	alfa = invtraincv * y0;
end