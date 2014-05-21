% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Apr 25 09:17:31 EDT 2014

% Define some covariance functions

TGDefineCovFuncs;

refyear=2010;
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

% asymptotic covariance associated with
% A exp((t0-t)/tau)
% so covariance of A^2 exp((2*t0-t1-t2)/tau)
%kASYMPT=@(years1,years2,thetas) thetas(1).^2*exp((2*thetas(2)-bsxfun(@plus,years1,years2'))/thetas(3));

%%%%%

modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 ...
                    thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2* ...
                                       (thetas(13)+thetas(14)*fp1fp2) ...
                                       + kGEOGG(ad,thetas([3 6 7]))) ...
    + kDP(t1,t2,1).*(kGEOGG(ad,thetas([8 9 10]))) + kDELTAG(ad,1).*(thetas(11).^2 + kDELTA(dt1t2,thetas(12)));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
100 1 1e4 % GIS amplitude
100 1 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 1.5 9.5 % temporal roughness (needs to be >= 1.5 for first derivative to be defined)

5 1 20 % geographic length scale
2.5 .5 9.5 % geographic roughness

.5 1e-3 100 % linear regional amplitude
5 .5 15 %  linear regional length scale
2.5 .5 9.5 % linear geographic roughness

100 1e-2 1e4 % local offset
10 1e-2 1e4 % white noise

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[length(modelspec(1).thet0)-[1 0]]; % fix GIS switches
modelspec(1).sublength=[6 9];

modelspec(1).subamp = [1 2 3 8 11 12];
modelspec(1).subamplinear = [8 11];
modelspec(1).subampglobal = [1];
modelspec(1).subamplocal = [11 12];
modelspec(1).subampnonlinear = [1 2 3 12];
modelspec(1).subampoffset = [11];
modelspec(1).subampGIS = [2];
modelspec(1).subampnoise = [12];
modelspec(1).label='default';

