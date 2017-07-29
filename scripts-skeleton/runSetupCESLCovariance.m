% set up covariance structures
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-07-25 18:13:30 -0400

clear modelspec;

% Define some covariance functions

refyear=2010;
CESLDefineCovFuncs;

% define some basic covariance function contributors
% using the different types of covariance functions defined
% in CESLDefineCovFuncs

cvfunc.G = @(dt1t2,thetas) kMat3(dt1t2,thetas(1:2));
cvfunc.L = @(t1,t2,ad,thetas) kDP(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
cvfunc.M = @(dt1t2,ad,thetas) kMat3(dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
cvfunc.I = @(dt1t2,ad,thetas,fp1fp2) kMat3(dt1t2,thetas(1:2)) .* fp1fp2;
cvfunc.W = @(dt1t2,ad,thetas) kDELTA(dt1t2,thetas(1)) .* kDELTAG(ad,1);
cvfunc.f0 = @(ad,thetas) kDELTAG(ad,thetas(1));

% first derivatives of the above

dcvfunc.G = @(t1,t2,dt1t2,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2));
dcvfunc.L = @(t1,t2,ad,thetas) kDPd(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
dcvfunc.M = @(t1,t2,dt1t2,ad,thetas) kMat3d(t1,t2,dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
dcvfunc.I = @(t1,t2,dt1t2,ad,thetas,fp1fp2) kMat3d(t1,t2,dt1t2,thetas(1:2)) .* fp1fp2;
dcvfunc.W = 0;
dcvfunc.f0 = 0;

% second derivatives of the above

ddcvfunc.G = @(dt1t2,thetas) kMat3dd(dt1t2,thetas(1:2));
ddcvfunc.L = @(t1,t2,ad,thetas) kDPdd(t1,t2,thetas(1)) .* kMat1(ad,[1 thetas(2)]) .* (ad<360);
ddcvfunc.M = @(dt1t2,ad,thetas) kMat3dd(dt1t2,thetas(1:2)) .* kMat1(ad,[1 thetas(3)]) .* (ad<360);
ddcvfunc.I = @(dt1t2,ad,thetas,fp1fp2) kMat3dd(dt1t2,thetas(1:2)) .* fp1fp2;
ddcvfunc.W = 0;
ddcvfunc.f0 = 0;

%%%%%

ii=1;

% 1. Global + Linear + Long-Wave Matern + Short-Wave Matern + local offset

modelspec(ii).label = 'GLMM';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.M(dt1t2,ad,thetas(8:10)) +  cvfunc.f0(ad,thetas(11));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(8:10));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.M(dt1t2,ad,thetas(8:10));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global amplitude
500 100 3e4 % time scale

1 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
500 100 3e3 % time scale
5 .5 20 % geographic length scale

100 1 1e4 % regional amplitude
500 10 3e3 % time scale
.1 .01 .5 % geographic length scale

100 1e-2 1e4 % local offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[4 7 10];

modelspec(ii).subamp = [1 3 5 8 11];
modelspec(ii).subamplinear = [3];
modelspec(ii).subampglobal = [1];
modelspec(ii).subampoffset = [11];
modelspec(ii).subampregmat = [5 8];
modelspec(ii).subampregmat1 = [5];
modelspec(ii).subampregmat2 = [8];
modelspec(ii).subampnoise = [];

modelspec(ii).subHPlinear = [3 4];
modelspec(ii).subHPglobal = [1 2];
modelspec(ii).subHPoffset = [11];
modelspec(ii).subHPregmat = [5 6 7 8 9 10];
modelspec(ii).subHPregmat1 = [5 6 7];
modelspec(ii).subHPregmat2 = [8 9 10];
modelspec(ii).subHPnoise = [];
modelspec(ii).subHPtsregmat = [6 9];
modelspec(ii).subHPtsglobal = [2];

%%%%%

%%%%

ii=ii+1;

% 1. Global + Linear + Long-Wave Matern + Short-Wave Matern (same time scale) + white noise + local offset


modelspec(ii).label = 'GLMmW';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.M(dt1t2,ad,thetas([8 6 9])) + cvfunc.W(dt1t2,ad,thetas(10)) + cvfunc.f0(ad,thetas(11));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 6 9]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.M(dt1t2,ad,thetas([8 6 9]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global amplitude
500 100 3e4 % time scale

1 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
500 100 3e3 % time scale
5 .5 20 % geographic length scale

100 1 1e4 % regional amplitude
.1 .01 .5 % geographic length scale

10 1e-2 1e4 % white noise
    
100 1e-2 1e4 % local offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[4 7];

modelspec(ii).subamp = [1 3 5 8 10 11];
modelspec(ii).subamplinear = [3];
modelspec(ii).subampglobal = [1];
modelspec(ii).subampoffset = [11];
modelspec(ii).subampregmat = [5 8];
modelspec(ii).subampregmat1 = [5];
modelspec(ii).subampregmat2 = [8];
modelspec(ii).subampnoise = [10];

modelspec(ii).subHPlinear = [3 4];
modelspec(ii).subHPglobal = [1 2];
modelspec(ii).subHPoffset = [11];
modelspec(ii).subHPregmat = [5 6 7 8 9];
modelspec(ii).subHPregmat1 = [5 6 7];
modelspec(ii).subHPregmat2 = [8 9];
modelspec(ii).subHPnoise = [10];
modelspec(ii).subHPtsregmat = [6];
modelspec(ii).subHPtsglobal = [2];


%%%%%

ii=ii+1;

% 3. Global Matern + Regional Linear + Regional Matern + White Noise (GLMW)

modelspec(ii).label = 'GLMW';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9)) +thetas(10)^2;

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global amplitude
500 100 3e4 % time scale

1 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
500 100 3e3 % time scale
5 .5 20 % geographic length scale

10 1e-2 1e4 % white noise
    
100 1e-2 1e4 % local offset
10 0 1e4 % global offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[4 7];

modelspec(ii).subamp = [1 3 5 8 9 10];
modelspec(ii).subamplinear = [3];
modelspec(ii).subampglobal = [1];
modelspec(ii).subampoffset = [9 10];
modelspec(ii).subampregmat = [5];
modelspec(ii).subampnoise = [8];

modelspec(ii).subHPlinear = [3 4];
modelspec(ii).subHPglobal = [1 2];
modelspec(ii).subHPoffset = [9 10];
modelspec(ii).subHPregmat = [5 6 7];
modelspec(ii).subHPnoise = [8];
modelspec(ii).subHPtsregmat = [6];
modelspec(ii).subHPtsglobal = [2];

%%%%%

ii = ii+1;

% Regional LF Matern + Regional MF Matern + Regional HF Matern +  White Noise (LMHW)

modelspec(ii).label = 'LMHW';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.M(dt1t2,ad,thetas(1:3)) + cvfunc.M(dt1t2,ad,thetas(4:6)) + cvfunc.M(dt1t2,ad,thetas([7 8 6])) + cvfunc.W(dt1t2,ad,thetas(9)) + cvfunc.f0(ad,thetas(10));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.M(t1,t2,dt1t2,ad,thetas(1:3)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(4:6)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([7 8 6]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.M(dt1t2,ad,thetas(1:3)) + ddcvfunc.M(dt1t2,ad,thetas(4:6)) + ddcvfunc.M(dt1t2,ad,thetas([7 8 6]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

5000 1 1e6 % regional amplitude
5000 2000 3e4 % time scale
5 .1 25 % geographic length scale

2000 1 1e6 % regional amplitude
500 100 2e3 % time scale
5 .1 25 % geographic length scale

2000 1 1e6 % regional amplitude
30 5 100 % time scale

100 1e-2 1e4 % white noise    
100 1e-2 1e4 % local offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[3 6];

modelspec(ii).subamp = [1 4 7 9 10];
modelspec(ii).subamplinear = [];
modelspec(ii).subampglobal = [];
modelspec(ii).subampoffset = [9];
modelspec(ii).subampregmat = [1 4 7];
modelspec(ii).subampnoise = [9];

modelspec(ii).subHPlinear = [];
modelspec(ii).subHPglobal = [];
modelspec(ii).subHPoffset = [10];
modelspec(ii).subHPregmat = [1:8];
modelspec(ii).subHPnoise = [9];
modelspec(ii).subHPtsregmat = [2 5 8];
modelspec(ii).subHPtsglobal = [];
