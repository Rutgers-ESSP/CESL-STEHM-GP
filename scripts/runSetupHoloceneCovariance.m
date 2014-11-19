% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Nov 18 11:47:02 EST 2014

% 1. GLMW
% 2. GrinstedGLMW - GLMW with global hyperparameters set to maximize likelihood of Grinsted curve
% 3. GLW
% 4. GMW
% 5. GW
% 6. LMW
% 7. LW
% 8. MW

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

% 1. Global Matern + Regional Linear + Regional Matern + GIS +  White Noise (GLMIW)

modelspec(1).label = 'GLMW';

modelspec(1).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9));

modelspec(1).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7));

modelspec(1).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

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

];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';

modelspec(1).subfixed=[];
modelspec(1).sublength=[4 7];

modelspec(1).subamp = [1 3 5 8 9];
modelspec(1).subamplinear = [3];
modelspec(1).subampglobal = [1];
modelspec(1).subampoffset = [9];
modelspec(1).subampregmat = [5];
modelspec(1).subampnoise = [8];

modelspec(1).subHPlinear = [3 4];
modelspec(1).subHPglobal = [1 2];
modelspec(1).subHPoffset = [9];
modelspec(1).subHPregmat = [5 6 7];
modelspec(1).subHPnoise = [8];


%%%%%%


%%%%%

% 3. Global Matern + Regional Linear + White Noise (GLW)

ii=3;
turnoff= [modelspec(1).subampregmat ];
freeze= [modelspec(1).subHPregmat ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GLW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%
% 4. GMW

ii=4;
turnoff= [modelspec(1).subamplinear ];
freeze= [modelspec(1).subHPlinear ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% 5. GW

ii=5;
turnoff= [modelspec(1).subampregmat modelspec(1).subamplinear ];
freeze= [modelspec(1).subHPregmat modelspec(1).subHPlinear ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='GW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);



% 6. Regional Linear + Regional Matern + White Noise (LMW)

ii=6;
turnoff= [modelspec(1).subampglobal ];
freeze= [modelspec(1).subHPglobal ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='LMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);



%%%%%

% 7. Regional Linear + White Noise (LW)

ii=7;
turnoff= [modelspec(1).subampglobal modelspec(1).subampregmat ];
freeze= [modelspec(1).subHPglobal modelspec(1).subHPregmat ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='LW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% 7. MW

ii=7;
turnoff= [modelspec(1).subampglobal modelspec(1).subamplinear ];
freeze= [modelspec(1).subHPglobal modelspec(1).subHPlinear ];

modelspec(ii) = modelspec(1);
modelspec(ii).label='MW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);


% 2. GrinstedGLMW

ii=2;
turnoff= [ ];
freeze= [modelspec(1).subHPglobal ];

ms = modelspec(4); % use GW
[thetGrin]= OptimizeHoloceneCovariance(Grinsted,ms,[2.4 2.0],[],[],.01);

modelspec(ii) = modelspec(1);
modelspec(ii).label='GrinstedGLMW';
%modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).thet0(modelspec(1).subHPglobal)=thetGrin(modelspec(1).subHPglobal);
