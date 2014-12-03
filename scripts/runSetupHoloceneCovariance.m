% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Dec 02 08:25:31 EST 2014

% 1. GLMW
% 2. GLMW-Grinsted - GLMW with global hyperparameters set to maximize likelihood of Grinsted curve
% 3. GLMW-NC
% 4. GLMW-1ts
% 5. GLMW-`1amp
% 5. GLW
% 6. GMW
% 7. GW
% 8. LMW
% 9. LW
% 10. MW

clear modelspec;

% Define some covariance functions

refyear=2010;
CESLDefineCovFuncs;

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

ii=1;

% 1. Global Matern + Regional Linear + Regional Matern + GIS +  White Noise (GLMIW)

modelspec(ii).label = 'GLMW';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9)) +thetas(10)^2;

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
1e3 1 3e4 % time scale

.5 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
1e3 1 3e3 % time scale
5 1 50 % geographic length scale

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


%%%%%%

% Grinsted prior

ii=2;
turnoff= [ ];
freeze= [modelspec(1).subHPglobal ];

clear ms;
ms.cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.W(dt1t2,ad,thetas(3));
ms.dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2));
ms.ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2));
ms.traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) ms.cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;
tluTGG=[
100 1 1e4 % global amplitude
200 1 3e4 % time scale
100 1e-3 1e4 % white noise
       ];
ms.thet0=tluTGG(:,1)';
ms.lb = tluTGG(:,2)';
ms.ub = tluTGG(:,3)';
ms.subfixed=[];
ms.sublength=[];
[thetGrin]= OptimizeHoloceneCovariance(Grinsted,ms,[2.4 2.0],[],[],.01);

modelspec(ii) = modelspec(1);
modelspec(ii).label='GLMW-GrinstedG';
%modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).thet0(modelspec(ii).subHPglobal)=thetGrin(1:2);


%%%%

% North Carolina prior

ii=ii+1;

wdataset=datasets{1};
sitesub=find(strcmpi('North Carolina-Sand Point',wdataset.sitenames));
sitesub=union(sitesub,find(strcmpi('North Carolina-Tump Point',wdataset.sitenames)));
datsub=find(ismember(wdataset.datid,wdataset.siteid(sitesub)));


clear ms;
ms.cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,[thetas(3:4)])+cvfunc.W(dt1t2,ad,thetas(5));
ms.dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,[thetas(3:4)]);
ms.ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,[thetas(3:4)]);
ms.traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) ms.cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;
tluTGG=[
100 1 1e4 % global amplitude
200 1 3e4 % time scale
.5 1e-3 100 % linear regional amplitude
5 .1 50 % linear length scale
100 1e-3 1e4 % white noise
       ];
ms.thet0=tluTGG(:,1)';
ms.lb = tluTGG(:,2)';
ms.ub = tluTGG(:,3)';
ms.subfixed=[];
ms.sublength=[4];
thetNC=OptimizeHoloceneCovariance(SubsetDataStructure(wdataset,datsub,sitesub),ms,[ 3.4 3.0],-1000,2000,.01);
modelspecNC=ms;

modelspec(ii) = modelspec(1);

turnoff= [ ];
freeze= [modelspec(1).subHPglobal ];

modelspec(ii).label='GLMW-NCG';
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).thet0(modelspec(ii).subHPglobal)=thetNC(1:2);

%%%%%


% Single timescale prior

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subHPtsregmat ];
freeze= [modelspec(ii).subHPtsregmat ];

modelspec(ii).label='GLMW-1ts';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas([5 2 7])) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9)) +thetas(10)^2;

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([5 2 7]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas([5 2 7]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

%%%%%


% Single timescale and amplitude prior

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampregmat modelspec(ii).subHPtsregmat ];
freeze= [modelspec(ii).subampregmat modelspec(ii).subHPtsregmat ];

modelspec(ii).label='GLMW-1amp1ts';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas([1 2 7])) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9))  + thetas(10)^2;;

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas([1 2 7]));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas([1 2 7]));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

%%%%


% 3. Global Matern + Regional Linear + White Noise (GLW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampregmat ];
freeze= [modelspec(ii).subHPregmat ];

modelspec(ii).label='GLW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%
% 4. GMW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPlinear ];

modelspec(ii).label='GMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% 5. GW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampregmat modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPregmat modelspec(ii).subHPlinear ];

modelspec(ii).label='GW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);



% 6. Regional Linear + Regional Matern + White Noise (LMW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal ];
freeze= [modelspec(ii).subHPglobal ];

modelspec(ii).label='LMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);



%%%%%

% 7. Regional Linear + White Noise (LW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal modelspec(ii).subampregmat ];
freeze= [modelspec(ii).subHPglobal modelspec(ii).subHPregmat ];

modelspec(ii).label='LW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% 7. MW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPglobal modelspec(ii).subHPlinear ];

modelspec(ii).label='MW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

