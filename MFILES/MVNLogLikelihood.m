function [logp,invcv] = MVNLogLikelihood(y,mu,Sigma,invcv)

% [logp,invcv] = MVNLogLikelihood(y,mu,Sigma,[invcv])
%
% Calculate log-likelihood for y from multivariate normal distribution with mean mu and covariance Sigma.
%
% mu should be a column vectors. Each sample in y should be a column.
%
% Output:
%
%   logp: log likelihood
%   invcv: Stored intermediary variables. Pass along to avoid inverting covariance matrix multiple times.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Jul 24 15:01:46 JST 2015

%%%%%
    
    errorflags=0;
    defval('doChol',1);
    defval('invcv',[]);
    y0=bsxfun(@minus,y,mu);


    if length(invcv)==0
        if doChol
            try
                L=chol(Sigma,'lower');
                %                    invSigma = L'\(L\eye(size(L)));
                alfa=L'\(L\y0);
            catch
                disp('Not positive definite!')
                errorflags=-1;
                [alfa,invSigma,s,r]=svdinv(Sigma,y0);
                doChol = 0;
            end
        else
            [alfa,invSigma,s,r]=svdinv(Sigma,y0);
        end
        
        if nargout>1
            if exist('L','var')
                invcv.L=L;
            end
            if exist('invSigma','var')
                invcv.invSigma=invSigma;
            end
            if exist('s','var')
                invcv.s=s;
            end
            if exist('r','var')
                invcv.r=r;
            end        
            invcv.alfa=alfa;
        end
    else
        if isfield(invcv,'L')
            doChol=1;
            L=invcv.L;
        else
            doChol=0;
            s=invcv.s;
            r=invcv.r;
            invSigma=invcv.invSigma;
        end
        alfa=invcv.alfa;
    end
    
    logptermsA = diag(-.5*abs(y0'*alfa))';
    logptermsB(2) = -.5*length(y0)*log(2*pi);
    % logptermsB(1) = -.5 * log(det(Sigma));
    if doChol
        logptermsB(1) = - sum(log(diag(L)));
    else
        logptermsB(1) = -.5 * sum(log((s(1:r))));
    end

    logp=sum(logptermsB)+logptermsA;

    %	logp = -.5*y0'*alfa - .5*log(prod(s(1:r)))- .5*length(y0)*log(2*pi);

end

function [alfa,invSigma,s,r]=svdinv(Sigma,y0)
    [m,n] = size(Sigma);
    [U,S,V] = svd(Sigma,0);
    s = diag(S);
    tol = max(m,n) * eps(max(s));
    r = sum(s > tol);
    invSigma = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
    alfa = invSigma * y0;
end