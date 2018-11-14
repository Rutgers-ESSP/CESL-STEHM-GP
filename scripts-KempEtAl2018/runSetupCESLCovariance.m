% set up covariance structures
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:08:37 -0400

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

clear modelspec;

ii=1;

modelspec(ii).label = 'GLMM-3ts';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.M(dt1t2,ad,thetas([8 9 10])) + cvfunc.W(dt1t2,ad,thetas(11)) + cvfunc.f0(ad,thetas(12));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 9 10]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.M(dt1t2,ad,thetas([8 9 10]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global amplitude
300 100 3e4 % time scale

1 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
300 100 3e4 % time scale
5 1 40 % geographic length scale

100 1 1e4 % regional amplitude
300 100 3e4 % time scale
.1 .01 1 % geographic length scale

10 1e-2 1e4 % white noise
    
100 1e-2 1e4 % local offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[4 7 10];

modelspec(ii).subamp = [1 3 5 8 11 12];
modelspec(ii).subamplinear = [3];
modelspec(ii).subampglobal = [1];
modelspec(ii).subampoffset = [12];
modelspec(ii).subampregmat = [5 8];
modelspec(ii).subampregmat1 = [5];
modelspec(ii).subampregmat2 = [8];
modelspec(ii).subampnoise = [11];

modelspec(ii).subHPlinear = [3 4];
modelspec(ii).subHPglobal = [1 2];
modelspec(ii).subHPoffset = [12];
modelspec(ii).subHPregmat = [5 6 7 8 9 10];
modelspec(ii).subHPregmat1 = [5 6 7];
modelspec(ii).subHPregmat2 = [8 9 10];
modelspec(ii).subHPnoise = [11];
modelspec(ii).subHPtsregmat = [6];
modelspec(ii).subHPtsglobal = [2];

ii=2;

modelspec(ii)=modelspec(1);
modelspec(ii).label = 'GLMm-2ts';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.M(dt1t2,ad,thetas([8 6 10])) + cvfunc.W(dt1t2,ad,thetas(11)) + cvfunc.f0(ad,thetas(12));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 6 10]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.M(dt1t2,ad,thetas([8 6 10]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

modelspec(ii).thet0(9) = 0;
modelspec(ii).subfixed=[9];

ii=3;

modelspec(ii)=modelspec(2);
modelspec(ii).label = 'LMm-1ts';

modelspec(ii).thet0([1 9])=0;
modelspec(ii).subfixed=[1 2 9];

ii=4;
modelspec(ii)=modelspec(1);
modelspec(ii).label = 'gLMm-1ts';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas([1 6])) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.M(dt1t2,ad,thetas([8 6 10])) + cvfunc.W(dt1t2,ad,thetas(11)) + cvfunc.f0(ad,thetas(12));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas([1 6])) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 6 10]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas([1 6])) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7)) + ddcvfunc.M(dt1t2,ad,thetas([8 6 10]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

modelspec(ii).thet0([2 9])=0;
modelspec(ii).subfixed=[2 9];

ii=5;
modelspec(ii)=modelspec(1);
modelspec(ii).label = 'GLM-2ts';

modelspec(ii).thet0([8])=0;
modelspec(ii).subfixed=[8 9];

ii=6;
modelspec(ii)=modelspec(1);
modelspec(ii).label = 'GLmm-1ts';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas([5 2 7])) + cvfunc.M(dt1t2,ad,thetas([8 2 10])) + cvfunc.W(dt1t2,ad,thetas(11)) + cvfunc.f0(ad,thetas(12));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([5 2 7])) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 2 10]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas([5 2 7])) + ddcvfunc.M(dt1t2,ad,thetas([8 2 10]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;


modelspec(ii).thet0([6 9])=0;
modelspec(ii).subfixed=[6 9];

ii=7;
modelspec(ii)=modelspec(6);
modelspec(ii).label = 'GLm-1ts';

modelspec(ii).thet0([6 8 9])=0;
modelspec(ii).subfixed=[6 8 9];

ii=8;

modelspec(ii).label = 'GMML-2tsx3ss';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.G(dt1t2,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas([5 2 6])) + cvfunc.M(dt1t2,ad,thetas([7 4 6])) + ...
    cvfunc.M(dt1t2,ad,thetas([8 2 9])) + cvfunc.M(dt1t2,ad,thetas([10 4 9])) + cvfunc.L(t1,t2,ad,thetas(11:12)) + cvfunc.W(dt1t2,ad,thetas(13)) + cvfunc.f0(ad,thetas(14));

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.G(t1,t2,dt1t2,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([5 2 6])) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([7 4 6])) + ...
    dcvfunc.M(t1,t2,dt1t2,ad,thetas([8 2 9])) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([10 4 9])) + dcvfunc.L(t1,t2,ad,thetas(11:12));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.G(dt1t2,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas([5 2 6])) + ddcvfunc.M(dt1t2,ad,thetas([7 4 6])) + ...
    ddcvfunc.M(dt1t2,ad,thetas([8 2 9])) + ddcvfunc.M(dt1t2,ad,thetas([10 4 9])) + ddcvfunc.L(t1,t2,ad,thetas(11:12));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;


tluTGG = [

50 1 1e4 % global amplitude 1
150 100 400 % time scale 1

100 1 1e4 % global amplitude 2
600 400 3e4 % time scale 2

50 1 1e4 % regional amplitude 1a
5 1 40 % geographic length scale 1
100 1 1e4 % regional amplitude 1b

50 1 1e4 % regional amplitude 2a
.5 .01 1 % geographic length scale 2
100 1 1e4 % regional amplitude 2b

1 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

10 1e-2 1e4 % white noise
    
100 1e-2 1e4 % local offset

];

modelspec(ii).thet0=tluTGG(:,1)';
modelspec(ii).lb = tluTGG(:,2)';
modelspec(ii).ub = tluTGG(:,3)';

modelspec(ii).subfixed=[];
modelspec(ii).sublength=[6 9 12];

modelspec(ii).subamp = [1 3 5 7 8 10 11 13 14];
modelspec(ii).subamplinear = [11];
modelspec(ii).subampglobal = [1 3];
modelspec(ii).subampoffset = [14];
modelspec(ii).subampregmat = [5 7 8 10];
modelspec(ii).subampregmat1 = [5 7];
modelspec(ii).subampregmat2 = [8 10];
modelspec(ii).subampnoise = [13];

modelspec(ii).subHPlinear = [11 12];
modelspec(ii).subHPglobal = [1 2 3 4];
modelspec(ii).subHPoffset = [14];
modelspec(ii).subHPregmat = [5 6 7 8 9 10];
modelspec(ii).subHPregmat1 = [5 6 7];
modelspec(ii).subHPregmat2 = [8 9 10];
modelspec(ii).subHPnoise = [14];
modelspec(ii).subHPtsregmat = [2 4];
modelspec(ii).subHPtsglobal = [2 4];

ii=9;

modelspec(ii)=modelspec(8);
modelspec(ii).label = 'GMML-2tsx3ss_lb';

modelspec(ii).lb([2 4])=[2 105];
modelspec(ii).ub(2)=105;

