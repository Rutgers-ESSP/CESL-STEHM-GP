% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Apr 23 08:20:08 EDT 2014



% Define some covariance functions

TGDefineCovFuncs;
defval('thetTG',[1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826]);


refyear=2010;
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

% asymptotic covariance associated with
% A exp((t0-t)/tau)
% so covariance of A^2 exp((2*t0-t1-t2)/tau)
%kASYMPT=@(years1,years2,thetas) thetas(1).^2*exp((2*thetas(2)-bsxfun(@plus,years1,years2'))/thetas(3));

%%%%%

%modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 thetas([5 6])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(12)+thetas(13)*fp1fp2) + kGEOG(ad,thetas([3 7])) + kDELTAG(ad,thetas(4))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([8 9])) + kDELTAG(ad,thetas(10))) + kDELTAG(ad,thetas(11));

modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(14)+thetas(15)*fp1fp2) + kGEOG(ad,thetas([3 6]))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([7 8])) + kDELTAG(ad,thetas(9))) + kDELTAG(ad,1).*(thetas(10).^2 + kMatG(dt1t2,thetas(11:13)));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
100 1 1e4 % GIS amplitude
100 1 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

5 1 20 % geographic length scale

.5 1e-2 100 % linear regional amplitude
5 .5 15 %  linear regional length scale
.5 1e-2 100 % linear local amplitude

100 1e-2 1e4 % local offset

10 1e-2 1e4 % local Matern
1e3 1 3e3
1.5 .5 9.5

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[length(modelspec(1).thet0)-[1 0]]; % fix GIS switches
modelspec(1).sublength=[6 8];

modelspec(1).subamp = [1 2 3 7 9 10 11];
modelspec(1).subamplinear = [7 9 10];
modelspec(1).subampglobal = [1 14];
modelspec(1).subamplocal = [9 10 11];
modelspec(1).subampnonlinear = [1 2 3 11];
modelspec(1).subampoffset = [10];
modelspec(1).subampGIS = [2];
modelspec(1).label='default';

%%%%%%%%

modelspec(2).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 ...
                    thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2* ...
                                       (thetas(15)+thetas(16)*fp1fp2) + kGEOG(ad,thetas([3 6]))) + kDP(t1,t2,1).*(thetas(7).^2 + kGEOG(ad,thetas([8 9])) + kDELTAG(ad,thetas(10))) + kDELTAG(ad,1).*(thetas(11).^2 + kMatG(dt1t2,thetas(12:14)));

modelspec(2).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(2).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
100 1 1e4 % GIS amplitude
100 1 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

5 1 20 % geographic length scale

0.1 1e-3 100 % linear global amplitude
.5 1e-3 100 % linear regional amplitude
5 .5 15 %  linear regional length scale
.5 1e-3 100 % linear local amplitude

100 1e-2 1e4 % local offset

10 1e-2 1e4 % local Matern
1e3 1 3e3
1.5 .5 9.5

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(2).thet0=tluTGG(:,1)';
modelspec(2).lb = tluTGG(:,2)';
modelspec(2).ub = tluTGG(:,3)';
modelspec(2).subfixed=[length(modelspec(2).thet0)-[1 0]]; % fix GIS switches
modelspec(2).sublength=[6 9];

modelspec(2).subamp = [1 2 3 7 8 10 11 12];
modelspec(2).subamplinear = [8 10 11 7];
modelspec(2).subampglobal = [1 7 15];
modelspec(2).subamplocal = [10 11 12];
modelspec(2).subampnonlinear = [1 2 3 12];
modelspec(2).subampoffset = [11];
modelspec(2).subampGIS = [2];
modelspec(2).label='globallin';


%%%%%%%

modelspec(3).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 ...
                    thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2* ...
                                       (thetas(14)+thetas(15)*fp1fp2) ...
                                       + kGEOGG(ad,thetas([3 6 7]))) ...
    + kDP(t1,t2,1).*(thetas(8).^2 + kGEOGG(ad,thetas([9 10 11]))) + kDELTAG(ad,1).*(thetas(12).^2 + kDELTA(dt1t2,thetas(13)));

modelspec(3).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(3).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
100 1 1e4 % GIS amplitude
100 1 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

5 1 20 % geographic length scale
2.5 .5 9.5 % geographic roughness

0.1 1e-3 100 % linear global amplitude
.5 1e-3 100 % linear regional amplitude
5 .5 15 %  linear regional length scale
2.5 .5 9.5 % linear roughness

100 1e-2 1e4 % local offset
10 1e-2 1e4 % white noise

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(3).thet0=tluTGG(:,1)';
modelspec(3).lb = tluTGG(:,2)';
modelspec(3).ub = tluTGG(:,3)';
modelspec(3).subfixed=[length(modelspec(3).thet0)-[1 0]]; % fix GIS switches
modelspec(3).sublength=[6 9];

modelspec(3).subamp = [1 2 3 8 9 12 13];
modelspec(3).subamplinear = [9 12 8];
modelspec(3).subampglobal = [1 8];
modelspec(3).subamplocal = [12 13];
modelspec(3).subampnonlinear = [1 2 3 13];
modelspec(3).subampoffset = [12];
modelspec(3).subampGIS = [2];
modelspec(3).label='globalroughness';


%%%%%%%

modelspec(4).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 ...
                    thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2* ...
                                       (thetas(13)+thetas(14)*fp1fp2) ...
                                       + kGEOGG(ad,thetas([3 6 7]))) ...
    + kDP(t1,t2,1).*(kGEOGG(ad,thetas([8 9 10]))) + kDELTAG(ad,1).*(thetas(11).^2 + kDELTA(dt1t2,thetas(12)));

modelspec(4).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(4).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

100 1 1e4 % global non-GIS amplitude
100 1 1e4 % GIS amplitude
100 1 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

5 1 20 % geographic length scale
2.5 .5 9.5 % geographic roughness

.5 1e-3 100 % linear regional amplitude
5 .5 15 %  linear regional length scale
2.5 .5 9.5 % linear roughness

100 1e-2 1e4 % local offset
10 1e-2 1e4 % white noise

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(4).thet0=tluTGG(:,1)';
modelspec(4).lb = tluTGG(:,2)';
modelspec(4).ub = tluTGG(:,3)';
modelspec(4).subfixed=[length(modelspec(4).thet0)-[1 0]]; % fix GIS switches
modelspec(4).sublength=[6 9];

modelspec(4).subamp = [1 2 3 8 11 12];
modelspec(4).subamplinear = [8 11];
modelspec(4).subampglobal = [1];
modelspec(4).subamplocal = [11 12];
modelspec(4).subampnonlinear = [1 2 3 12];
modelspec(4).subampoffset = [11];
modelspec(4).subampGIS = [2];
modelspec(4).label='globalroughness2';

