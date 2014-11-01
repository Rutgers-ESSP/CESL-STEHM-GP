function [thet,logp,acceptedcount]=SLNIGPSamplePT(x0,y0,dx0,dy0,traincvaug,modelspec,thet0,lb,ub,spacex,temps,Nsamples,Nthin,Nburnin,savefile,plotfile)

% Parallel tempering MCMC
% note does not yet implement thinning/burnin
% currently using log step sizes (I think equivalent to a log-uniform prior)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 01 09:56:15 EDT 2014

    defval('spacex',ones(size(x0)));
    defval('Nsamples',1000);
    defval('Nthin',10);
    defval('Nburnin',10);
    defval('savefile','');
    defval('temps',1:4);
    defval('plotfile','');

    %stepsize=(log(ub)-log(lb))/1000;
    stepsize=ones(size(thet0))*0.01; % 1% change in value
    Ntemps=length(temps);
    
    testfunc = @logprobNI;
    
    thet=repmat(thet0,1,1,Ntemps);    
    logp(1,:) = repmat(testfunc(x0,y0,dx0,dy0,traincvaug,modelspec,thet(1,:,1),spacex),1,Ntemps);

    tic
    
    for nn=2:Nsamples
        thet(nn,:,:) = thet(nn-1,:,:);
        logpold=logp(nn-1,:);
        acceptedcount(nn,:)=zeros(1,Ntemps);
 
        % normal Metropolis step
        parfor pp=1:Ntemps
            wthet = thet(nn,:,pp);
            wlogpold = logpold(pp);
            for ii=randperm(length(thet0))       
                proposal = wthet;
                proposal(ii) = exp(log(wthet(ii)) + randn*stepsize(ii));
                logpproposed=testfunc(x0,y0,dx0,dy0,traincvaug,modelspec,proposal,spacex);
                if (log(rand) < ((logpproposed-logpold)*temps(pp))).*(sum(proposal<lb)==0).*(sum(proposal>ub)==0)
                    disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- accepted -- %0.3f (%0.3f)',[pp nn ii logpproposed wlogpold]));
                    wlogpold = logpproposed;
                    wthet = proposal;
                    acceptedcount(nn,pp) = acceptedcount(nn,pp)+1/length(thet0);
                else
                    disp(sprintf('%0.0f -- Round %0.0f - %0.0f -- rejected -- %0.3f (%0.3f)',[pp nn ii logpproposed wlogpold]));
                end
            end
            logp(nn,pp)=wlogpold;
            thet(nn,:,pp) = wthet;
        end
        
        % now mix
        mixers=randperm(Ntemps);
        for ppp=1:2:length(mixers)
            m1=mixers(ppp); m2=mixers(ppp+1);
            logpmix = temps(m1)*logp(nn,m2) + temps(m2)*logp(nn,m1) - temps(m1)*logp(nn,m1) - temps(m2)*logp(nn,m2);
            if log(rand) < logpmix
                mixed(nn,[0 1]+ppp) = [m1 m2];
                thetm1 = thet(nn,:,m1);
                thet(nn,:,m1)=thet(nn,:,m2);
                thet(nn,:,m2) = thetm1;
                logpm1 = logp(nn,m1);
                logp(nn,m1)=logp(nn,m2);
                logp(nn,m2)=logpm1;
                disp(sprintf('Mixed chains %0.0f and %0.0f',[m1 m2]));
            else
                mixed(nn,[0 1]+ppp)=[0 0];
            end
        end
        
        if length(savefile)>0
            save(savefile,'thet','logp','acceptedcount','mixed');
        end
        if length(plotfile)>0
            if mod(nn,10)==0
                clf;
                plot(1:nn,logp);
                ylabel('log p');
                pdfwrite(plotfile);
            end
            
        end
        timepassed=toc;
        disp(sprintf('%0.0f seconds passed -- %0.0f estimated seconds remaining',[timepassed timepassed/(nn-1)*(Nsamples-nn)]));
        
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
