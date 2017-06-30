% set up covariance structures
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-06-30 15:26:23 -0400

% Alternative covariance structures:
%
% 1. GLMW [known as ML22 in paper]
% 2. GLMW-1ts [known as ML21 in paper]
% 3. GLMW-1amp1ts [known as ML11 in paper]
% 4. GLW
% 5. GMW
% 6. GW
% 7. LMW
% 8. LW
% 9. MW

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

% 1. Global Matern + Regional Linear + Regional Matern +  White Noise (GLMW)

modelspec(ii).label = 'GLMW';

modelspec(ii).cvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) cvfunc.G(dt1t2,thetas(1:2)) + cvfunc.L(t1,t2,ad,thetas(3:4)) + cvfunc.M(dt1t2,ad,thetas(5:7)) + cvfunc.W(dt1t2,ad,thetas(8)) + cvfunc.f0(ad,thetas(9)) +thetas(10)^2;

modelspec(ii).dcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) dcvfunc.G(t1,t2,dt1t2,thetas(1:2)) + dcvfunc.L(t1,t2,ad,thetas(3:4)) + dcvfunc.M(t1,t2,dt1t2,ad,thetas(5:7));

modelspec(ii).ddcvfunc =  @(t1,t2,dt1t2,thetas,ad,fp1fp2) ddcvfunc.G(dt1t2,thetas(1:2)) + ddcvfunc.L(t1,t2,ad,thetas(3:4)) + ddcvfunc.M(dt1t2,ad,thetas(5:7));

modelspec(ii).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(ii).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
500 100 3e4 % time scale

.5 1e-3 100 % linear regional amplitude
5 .5 50 %  linear regional length scale

100 1 1e4 % regional amplitude
500 100 3e3 % time scale
5 1 20 % geographic length scale

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


% Global Matern + Regional Linear + White Noise (GLW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampregmat ];
freeze= [modelspec(ii).subHPregmat ];

modelspec(ii).label='GLW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

%%
% GMW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPlinear ];

modelspec(ii).label='GMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% GW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampregmat modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPregmat modelspec(ii).subHPlinear ];

modelspec(ii).label='GW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);



% Regional Linear + Regional Matern + White Noise (LMW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal ];
freeze= [modelspec(ii).subHPglobal ];

modelspec(ii).label='LMW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

modelspec(ii).lb(modelspec(ii).subampregmat) = 1; % reduce lower bound so can assess where the optimum is with minimal constraint

%%%%%

% Regional Linear + White Noise (LW)

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal modelspec(ii).subampregmat ];
freeze= [modelspec(ii).subHPglobal modelspec(ii).subHPregmat ];

modelspec(ii).label='LW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

% MW

ii=ii+1;
modelspec(ii) = modelspec(1);

turnoff= [modelspec(ii).subampglobal modelspec(ii).subamplinear ];
freeze= [modelspec(ii).subHPglobal modelspec(ii).subHPlinear ];

modelspec(ii).label='MW';
modelspec(ii).thet0(turnoff)=0;
modelspec(ii).subfixed=union(modelspec(ii).subfixed,freeze);

