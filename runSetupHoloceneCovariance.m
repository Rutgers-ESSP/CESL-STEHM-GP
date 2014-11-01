% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Oct 27 13:33:05 EDT 2014

% Define some covariance functions

TGDefineCovFuncs;

refyear=2010;
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

% kMat3 = @(dx,thetas)thetas(1).^2.*(1+sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
% kMat1 = @(dx,thetas)thetas(1).^2.*(1).*exp(-dx/thetas(2))

kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1);
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1));


cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.L = @(t1,t2,ad,thetas) kDP(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
cvfunc.M = @(dt1t2,ad,thetas) kMat3(dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
cvfunc.I = @(dt1t2,ad,thetas,fp1fp2) kMat3(dt1t2,thetas(1:2)) .* fp1fp2;
cvfunc.W = @(dt1t2,ad,thetas) kDELTA(dt1t2,thetas(1)) .* kDELTAG(ad,1);
cvfunc.f0 = @(ad,thetas) kDELTAG(ad,thetas(1));

dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.L = @(t1,t2,ad,thetas) kDPd(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
dcvfunc.M = @(t1,t2,dt1t2,ad,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
dcvfunc.I = @(t1,t2,dt1t2,ad,thetas,fp1fp2) kMat3d(t1,t2,dt1t2,thetas(1:2)) .* fp1fp2;
dcvfunc.W = 0;
dcvfunc.f0 = 0;

ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.L = @(t1,t2,ad,thetas) kDPdd(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
ddcvfunc.M = @(dt1t2,ad,thetas) kMat3dd(dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
ddcvfunc.I = @(dt1t2,ad,thetas,fp1fp2) kMat3dd(dt1t2,thetas(1:2)) .* fp1fp2;
ddcvfunc.W = 0;
ddcvfunc.f0 = 0;

%%%%%

% 1. Global Matern + Regional Linear + Regional Matern + GIS +  White Noise (GLMIW)

modelspec(1).label = 'GLMIW';

modelspec(1).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.I(dt1t2,ad,thetas(8:9),fp1fp2) + cvfunc.W(dt1t2,ad,thetas(10)) + cvfunc.f0(ad,thetas(11));

modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.I(t1,t2,dt1t2,ad,thetas(8:9),fp1fp2);

modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.I(dt1t2,ad,thetas(8:9),fp1fp2);

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
1e3 1 3e3 % time scale

.5 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
1e3 1 3e3 % time scale
5 1 50 % geographic length scale

100 1 1e4 % GIS amplitude
1e3 1 3e3 % time scale

10 1e-2 1e4 % white noise
    
100 1e-2 1e4 % local offset

];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';

modelspec(1).subfixed=[];
modelspec(1).sublength=[4 7];

modelspec(1).subamp = [1 3 5 8 10 11];
modelspec(1).subamplinear = [3];
modelspec(1).subampglobal = [1];
modelspec(1).subampoffset = [11];
modelspec(1).subampGIS = [8];
modelspec(1).subampregmat = [5];
modelspec(1).subampnoise = [10];

modelspec(1).subHPlinear = [3 4];
modelspec(1).subHPglobal = [1 2];
modelspec(1).subHPoffset = [11];
modelspec(1).subHPGIS = [8 9];
modelspec(1).subHPregmat = [5 6 7];
modelspec(1).subHPnoise = [10];

%%%%%%

% 2. Regional Linear + Regional Matern + White Noise (LMW)

ii=2;
turnoff= [modelspec(1).subampglobal modelspec(1).subampGIS];
freeze= [modelspec(1).subHPglobal modelspec(1).subHPGIS];

modelspec(ii) = modelspec(1);
modelspec(ii).label='LMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%%%%%

% 3. Regional Linear + Regional Matern + GIS + White Noise (LMIW)

ii=3;
turnoff= [modelspec(1).subampglobal];
freeze= [modelspec(1).subHPglobal];

modelspec(ii) = modelspec(1);
modelspec(ii).label='LMIW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%%%%

% 4. Global Matern + Regional Linear + White Noise (GLW)

ii=4;
turnoff= [modelspec(1).subampregmat modelspec(1).subampGIS];
freeze= [modelspec(1).subHPregmat modelspec(1).subHPGIS];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GLW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%%%%%%%%

% 5. Global Matern + Regional Matern + Regional Linear + White Noise (GLMW)

ii=5;
turnoff= [modelspec(1).subampGIS];
freeze= [modelspec(1).subHPGIS];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GLMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%%%%

% 6. Global Matern + GIS + Regional Linear + White Noise (GLIW)

ii=6;
turnoff= [modelspec(1).subampregmat];
freeze= [modelspec(1).subHPregmat];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GLIW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);
