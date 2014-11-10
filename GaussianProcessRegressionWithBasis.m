function [basisest,basisestcv,f,V,fbasis,Vbasis,ferr,Verr,logp,alfa,errorflags,invtraincv,invcv] = GaussianProcessRegressionWithBasis(x0,y0,x,traincv,cvfunc,testcv2,basisfunc,basistest,basismean,basiscv,invcv)

% [basisest,basisestcv,f,V,fbasis,Vbasis,ferr,Verr,logp,alfa,errorflags,invtraincv,invcv] = GaussianProcessRegressionWithBasis(x0,y0,x,traincv,cvfunc,[testcv2],[basisfunc],[basistest],[basismean],[basiscv],[invcv])
%
% Modeled on Rasmussen & Williams (2006)
%
% INPUT
%
% 	x0			x at training points
%	y0			y at training points (demeaned)
%	x			x at which to approximate
%	traincv		covariance matrix among covariance points
%	cvfunc		EITHER a handle for a function that takes x and x0 as inputs
%				OR the covariance matrix between x and x0
%	testcv2		IF cvfunc is passed as a matrix, then the covariance matrix among x
%   basisfunc   EITHER a handle for a function that takes x as inputs
%               and returns explicit basis functions (each function in a column)
%               OR the values of the explicit basis function at each training x
%   basistest   IF basisfunc is passed as a matrix, then the explicit basis function
%               at each test point x
%   basismean   mean of the Gaussian prior for the coefficients of the explicit basis functions
%               (column vector) (set = [] for vague prior)
%   basiscv     prior covariance for the coefficients of the explicit basis functions (if ndim = 1, interpreted as standard deviations)
%
% Assumes a zero-mean Gaussian process; add/subtract the means separately otherwise.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 09 17:05:10 EST 2014
    
%%%%%
    
    errorflags=0;
    defval('tol',1e-6);
    defval('doChol',0);
    defval('invcv',[]);
    defval('testcv2',[]);
    defval('basisfunc',[]);
    defval('basistest',[]);
    defval('basismean',[]);
    defval('basiscv',[]);

    dobasis = length(basisfunc)>0;
    dovagueprior = length(basismean)==0;
    reportinvcv = (nargout>8);
    reportinvtraincv = (nargout>7);
    dologp = (nargout>4);
    
    if length(basismean)>1
        basismean=basismean(:);
    end
    if prod(size(basiscv)) == length(basiscv)
        basiscv = diag(basiscv.^2);
    end
    
    
    
    [f0,V0,logp,alfa,errorflags,invtraincv,invcv] = GaussianProcessRegression(x0,y0,x,traincv,cvfunc,testcv2,invcv);
    if length(testcv2)==0
        testcv = feval(cvfunc,x0,x);
    else
        testcv = cvfunc;
    end   
    invtraincvtestcv = invtraincv * testcv;
 
    % basisest = (Binv + basistrain * invtraincv * basistrain') \ (basistrain * invtraincv * y0 + Binv * basismean);
    % R = basisest - basistrain * invtraincv * testcv;
    % f = f + R' * basisest
    % V = V + R' * ((Binv + basistrain * invtraincv * basistrain') \ R);
    %
    % for vague priors: Binv = 0
    
    
    if dobasis
        if length(basistest) == 0
            basistrain = feval(basisfunc,x0)';
            basistest = feval(basisfunc,x)';
        else
            basistrain = basisfunc;
        end
        
        if dovagueprior
            basismean = zeros(size(basistrain,1),1);
        end
        
        
        calcBinv = 0;
        if length(invcv)==0; calcBinv=1; elseif ~isfield(invcv,'Binv'); calcBinv=1; else; Binv=invcv.Binv; end
        L2=[];
    
        if calcBinv
            if dovagueprior
                Binv=zeros(size(basistrain,1));
            else               
                if doChol
                    try
                        L2=chol(basiscv,'lower');
                        Binv = L2'\(L2\eye(size(L2)));
                    catch
                        disp('Basis prior CV not positive definite!')
                        errorflags=errorflags-.1;
                        [~,Binv,Bs,Br]=svdinv(basiscv,1);
                    end
                else
                    [~,Binv,Bs,Br]=svdinv(basiscv,1);
                end
            end

            if reportinvcv
                invcv.Binv = Binv;
            end
        end
        
        basisprec = (Binv + basistrain * invtraincv * basistrain');
        basisest = basisprec \ (basistrain * alfa + Binv * basismean);
 
        Rremove = -basistrain*invtraincvtestcv;
        Rbasis = basistest;
        R = Rbasis + Rremove;

        fadd = R' * basisest;
        Vadd =  R' * (basisprec \ R);
        f = f0 + fadd;
        V = V0 + Vadd;

        fadderr = Rremove' * basisest;
        Vadderr =  Rremove' * (basisprec \ Rremove);
        ferr = f0 + fadderr;
        Verr = V0 + Vadderr;
        
        fbasis = Rbasis' * basisest;
        Vbasis= Rbasis' * (basisprec \ Rbasis);
        
        if (dologp)
            C = invtraincv * basistrain' * (basisprec \ basistrain) * invtraincv;
            logptermsaddl(1) = -.5*abs(y0'*C*y0);
            if ~dovagueprior
                if length(L2)==0
                    logptermsadd(2) = -sum(log(diag(L2)));
                else
                    logptermsadd(2) = -.5*sum(log((Bs(1:Br))));
                end
            end
            try
                L3 = chol(basisprec,'lower');
                logptermsadd(3) = -sum(log(diag(L3)));
            catch
                [m,n] = size(basisprec);
                [U3,S3,V3] = svd(basisprec,0);
                s3 = diag(S3);
                tol = max(m,n) * eps(max(s3));
                r3 = sum(s3 > tol);
                logpterms(3) = -.5*sum(log((s3(1:r3))));
            end
        end
        basisestcv = basistest'\((Vbasis)/basistest);

    else
        basisest=[];
        basisestcv=[];
    end
    
     
    V=.5*(V+V');
    if dologp; logp=logp+sum(logptermsadd); end

end

function [alfa,invtraincv,s,r]=svdinv(traincv,y0)
    [m,n] = size(traincv);
    [U,S,V] = svd(traincv,0);
    s = diag(S);
    tol = max(m,n) * eps(max(s));
    r = sum(s > tol);
    invtraincv = V(:,1:r)*diag(s(1:r).^-1)*U(:,1:r)';      
    alfa = invtraincv * y0;
end