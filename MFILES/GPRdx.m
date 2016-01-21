function [dK,df,d2f,yoffset,f0,dV] = GPRdx(x0,y0,dx0,dy0,cvfunc0,Nderivs,spacex,varargin)

% [dK,df,d2f,yoffset,f0,dV] = GPRdx(x0,y0,dx0,dy0,cvfunc,[Nderivs],[spacex])
% [dK,df,d2f,yoffset,f0,dV] = GPRdx(x0,y0,dx0,dy0,modelspec)
% 
% Calculates increment of training covariance matrix for noisy GP regression.
%
% Note that cvfunc should be specified has cvfunc(x1,x2) or cvfunc(x1,x2,r1,r2)
% -- this is a different format than used in some other places in the code.
% OptimizeHoloceneCovariance and RegressHoloceneDataSets handle this translation,
% but if you access this function directly, please be aware.
%
% OUTPUTS
%
%     dK: term to add to diagonal = (df^2 .* dx0.^2) or (diag(df) * dx0 * diag(df))
%     df: first derivatives
%     d2f: second derivatives
%     yoffset: offset of y due to d2f
%     f0: noise-free mean projection at x0 
%     dV: uncertainty on derivatives (if calculated analytically)
%
%    
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Dec 11 11:50:48 EST 2014

%%%%%

    defval('Nprompt',200);
    defval('Nderivs',1);
    defval('testcv',[]);
    defval('doParallel',1);
    defval('spacex',[]);
    if isstruct(cvfunc0)
        doanalytical=1;
    else
        doanalytical=0;
    end
    
    if doanalytical
        if length(spacex)~=length(x0)
            cvfunc = @(x1,x2,r1,r2) cvfunc0.cvfunc(x1,x2);
            dcvfunc = @(x1,x2,r1,r2) cvfunc0.dcvfunc(x1,x2);
            ddcvfunc = @(x1,x2,r1,r2) cvfunc0.ddcvfunc(x1,x2);
            spacex = ones(size(x0));
        else
            cvfunc = @(x1,x2,r1,r2) cvfunc0.cvfunc(x1,x2,r1,r2);
            dcvfunc = @(x1,x2,r1,r2) cvfunc0.dcvfunc(x1,x2,r1,r2);
            ddcvfunc = @(x1,x2,r1,r2) cvfunc0.ddcvfunc(x1,x2,r1,r2);
        end
        
        testcv = cvfunc(x0,x0,spacex,spacex);
        dtestcv = dcvfunc(x0,x0,spacex,spacex);
        ddtestcv = ddcvfunc(x0,x0,spacex,spacex);
        if min(size(dy0))==1
            traincv = testcv+diag(dy0).^2;
        else
            traincv = testcv + dy0;
        end
  
        %[f0,~,~,~,~,~,invcv] = GaussianProcessRegression(x0,y0,x0,traincv,testcv',testcv);
        f0=[];
        [df,dV]= GaussianProcessRegression(x0,y0,x0,traincv,dtestcv',ddtestcv);
        d2f=[];
        Ndim=1;
        yoffset=zeros(size(y0));
        Nderivs=1;
    else
        
        
        
        if length(spacex)~=length(x0)
            cvfunc = @(x1,x2,r1,r2) cvfunc0(x1,x2);
            spacex = ones(size(x0));
        else
            cvfunc = @(x1,x2,r1,r2) cvfunc0(x1,x2,r1,r2);
        end
        

        for i=1:length(varargin)
            if strcmpi(varargin{i},'Nprompt')
                Nprompt=varargin{i+1};
                i=i+1;
            end
            if strcmpi(varargin{i},'noparallel')
                doParallel=0;
            end
        end

        yoffset=zeros(size(y0));

        if min(size(dx0))==1
            Ndim = 1;
        else
            Ndim = 2;
            dx02D = dx0;
            dx0 = sqrt(diag(dx0));
        end

        if length(testcv)==0
            testcv = cvfunc(x0,x0,spacex,spacex);
        end
        if min(size(dy0))==1
            traincv = testcv+diag(dy0).^2;
        else
            traincv = testcv + dy0;
        end

        [f0,~,~,~,~,~,invcv] = GaussianProcessRegression(x0,y0,x0,traincv,testcv',testcv);
        df=zeros(length(y0),1);
        dK=zeros(length(y0),1);

        if doParallel
            parfor i=1:length(dx0)
                if mod(i,Nprompt)==0
                    disp(i);
                end
                if Nderivs==1
                    [df(i)] = CalcDerivative(x0,y0,dx0,cvfunc,i,invcv,spacex);
                else
                    [df(i),d2f(i)] = CalcDerivative(x0,y0,dx0,cvfunc,i,invcv,spacex);
                end
            end
        else
            for i=1:length(dx0)
                if mod(i,Nprompt)==0
                    disp(i);
                end
                if Nderivs==1
                    [df(i)] = CalcDerivative(x0,y0,dx0,cvfunc,i,invcv,spacex);
                else
                    [df(i),d2f(i)] = CalcDerivative(x0,y0,dx0,cvfunc,i,invcv,spacex);
                end
            end
        end
        
        df=df(:);
        if Nderivs>1
            d2f=d2f(:);
        else
            
            d2f=[];
        end
    end
    
    if Ndim==1
        dK=df.^2.*dx0.^2;
        if Nderivs>1
            dK = dK + .25*d2f.^2 .*dx0.^4;
            yoffset = .5*dx0.^2 .* d2f;
        end
    else
        dK = diag(df) * dx02D * diag(df);
        if Nderivs>1
            dK = dK + .25*diag(d2f)*dx02D*dx02D'*diag(d2f);
            yoffset = .5*dx0.^2 .* d2f;
        end
    end
end

%%%%

function [df,d2f] = CalcDerivative(x0,y0,dx0,cvfunc,i,invcv,spacex)
    
% have to inefficienntly calculate full testcv and testcv2 so as
% not to require recalculating of fingerprints etc.
    
    defval('invcv',[]);
    Nderivs = nargout;
    traincv=[];
    stepsz = min(.01,max(1e-3*x0(i),.01*dx0(i)));
    if dx0(i)>0
        deltax = sparse(length(x0),1);
        deltax(i) = stepsz;
        
        mx0 = x0+deltax;
        %        testcv = cvfunc(x0,mx0,spacex{:}); testcv2 = cvfunc(mx0,mx0,spacex{:});
        testcv = cvfunc(x0,mx0(i),spacex,spacex(i,:))'; testcv2=cvfunc(mx0(i),mx0(i),spacex(i,:),spacex(i,:));
        f1 = GaussianProcessRegression(x0,y0,mx0,traincv,testcv,testcv2,invcv);

        mx0 = x0-deltax;
        testcv = cvfunc(x0,mx0(i),spacex,spacex(i,:))'; testcv2=cvfunc(mx0(i),mx0(i),spacex(i,:),spacex(i,:));
        f2 = GaussianProcessRegression(x0,y0,mx0,traincv,testcv,testcv2,invcv);

        df = (f1-f2)/(2*stepsz);
        if Nderivs==1
            d2f=[];
        else
            mx0 = x0+2*deltax;
            testcv = cvfunc(x0,mx0(i),spacex,spacex(i,:))'; testcv2=cvfunc(mx0(i),mx0(i),spacex(i,:),spacex(i,:));
             f3 = GaussianProcessRegression(x0,y0,mx0,traincv,testcv,testcv2,invcv);

            mx0 = x0-2*deltax;
            testcv = cvfunc(x0,mx0(i),spacex,spacex(i,:))'; testcv2=cvfunc(mx0(i),mx0(i),spacex(i,:),spacex(i,:));
             f4 = GaussianProcessRegression(x0,y0,mx0,traincv,testcv,testcv2,invcv);

           df3 = (f3-f1)/stepsz;
            df4 = (f2-f4)/stepsz;
            d2f = (df3-df4)/(3*stepsz);
        end
    else
        df=eps; d2f=eps;
    end
end
