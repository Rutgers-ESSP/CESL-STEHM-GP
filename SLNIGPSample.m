function [thet,logp,acceptedcount]=SLNIGPSample(x0,y0,dx0,dy0,traincvaug,cvfunc,thet0,lb,ub,spacex,Nsamples,Nthin,Nburnin,savefile)

% note does not yet implement thinning/burnin
% currently using log step sizes (I think equivalent to a log-uniform prior)
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Oct 27 21:27:52 EDT 2014

    defval('thet0',[20 50 10 2 2 2 2 2 0 2 0]);
    defval('lb',[0 5 0 1e-3 0 0 0 0 0 1e-3 0]);
    defval('ub',[Inf Inf Inf 10 Inf Inf Inf Inf Inf Inf Inf]);
    defval('basisX',[]);
    defval('spacex',ones(size(x0)));
    defval('Nsamples',1000);
    defval('Nthin',10);
    defval('Nburnin',10);
    defval('savefile','');

    lb=max(1e-12,lb);
    ub=max(1e-12,ub);
    stepsize=(log(ub)-log(lb))/1000;
    
    testfunc = @logprobNI;
    
    thet(1,:) = thet0;
    logp(1) = testfunc(x0,y0,dx0,dy0,traincvaug,cvfunc,thet(1,:),spacex);
    acceptedcount(1)=1;
    for nn=2:Nsamples
        thet(nn,:) = thet(nn-1,:);
        logpold=logp(nn-1);
        acceptedcount(nn)=0;
        for ii=1:length(thet0)          
            proposal = thet(nn,:);
            proposal(ii) = exp(log(thet(nn,ii)) + randn*stepsize(ii));
            logpproposed=testfunc(x0,y0,dx0,dy0,traincvaug,cvfunc,proposal,spacex);
            acceptance = exp(logpproposed-logpold);
            if rand < acceptance
                logpold = logpproposed;
                thet(nn,:) = proposal;
                acceptedcount(nn) = acceptedcount(nn)+1/length(thet0);
                disp(sprintf('Round %0.0f - %0.0f -- accepted -- %0.3f',[nn ii logpproposed]));
            else
                disp(sprintf('Round %0.0f - %0.0f -- rejected -- %0.3f',[nn ii logpproposed]));
            end
        end
        logp(nn) = logpold;
        if length(savefile)>0
            save(savefile,'thet','logp','acceptedcount');
        end
    end
    
    
end

function logp=logprobNI(x0,y0,dx0,dy0,traincvaug,modelspec,thet,spacex)

    if isstruct(modelspec)
        mspec.cvfunc = @(x1,x2,r1,r2) modelspec.cvfunc(x1,x2,thet,r1,r2);
        mspec.dcvfunc = @(x1,x2,r1,r2) modelspec.dcvfunc(x1,x2,thet,r1,r2);
        mspec.ddcvfunc = @(x1,x2,r1,r2) modelspec.ddcvfunc(x1,x2,thet,r1,r2);
        cvfunc=modelspec.cvfunc;
    else
        mspec = @(x1,x2,r1,r2) modelspec(x1,x2,thet,r1,r2); % this is cvfunc
        cvfunc = modelspec;
    end
    
    if min(size(dy0))==1
        if sum(abs(dx0))>0
            [dK,df,d2f,yoffset] = GPRdx(x0,y0,dx0,diag(dy0.^2)+traincvaug(thet),mspec,2,spacex);
        else
            dK = 0;
        end
        
        logp = logprob(y0-yoffset,@(theta) cvfunc(x0,x0,theta,spacex,spacex)+diag(dy0.^2)+diag(dK)+traincvaug(theta),thet,[]);
    else
        if sum(abs(dx0))>0
            [dK,df,d2f,yoffset] = GPRdx(x0,y0,dx0,dy0+traincvaug(thet),mspec,2,spacex);
        else
            dK = 0;
        end
 
        logp = logprob(y0-yoffset,@(theta) cvfunc(x0,x0,theta,spacex,spacex)+dy0+diag(dK)+traincvaug(theta),thet,[]);
    end

end


function logp=logprob(y0,traincv,x,basisX)

    errorflags=0;
    
    tcv=traincv(x);
    try
        L=chol((tcv),'lower');
        %	catch
        %			disp('fall back');
        %			L=cholcov(tcv,'lower');
        %			keyboard
        %		end
        alfa=L'\(L\y0);
        doChol=1;
    catch
        disp('Not positive definite!')
        errorflags=-1;
        doChol = 0;
        [m,n] = size(tcv);
        [U,S,V] = svd(tcv,0);
        s = diag(S);
        tol = max(m,n) * eps(max(s));
        r = sum(s > tol);
        invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
        alfa = invtraincv * y0;
        %	            disp('SVD')
    end

    logpterms(1) = -.5*abs(y0'*alfa);
    logpterms(3) = -.5*length(y0)*log(2*pi);
    % logpterms(2) = -.5 * log(det(traincv));
    if doChol
        logpterms(2) = - sum(log(diag(L)));
    else
        logpterms(2) = -.5 * sum(log((s(1:r))));
    end
    logp=sum(logpterms);
    if errorflags == -1
        logp = -1e20;
    end
    
end
