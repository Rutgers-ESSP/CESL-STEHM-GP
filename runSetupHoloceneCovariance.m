% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Mar 4 11:12:37 EST 2014


% define some covariance functions

TGDefineCovFuncs;
defval('thetTG',[1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826]);


refyear=2010;
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));

% asymptotic covariance associated with
% A exp((t0-t)/tau)
% so covariance of A^2 exp((2*t0-t1-t2)/tau)
%kASYMPT=@(years1,years2,thetas) thetas(1).^2*exp((2*thetas(2)-bsxfun(@plus,years1,years2'))/thetas(3));


%modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 thetas([5 6])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(12)+thetas(13)*fp1fp2) + kGEOG(ad,thetas([3 7])) + kDELTAG(ad,thetas(4))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([8 9])) + kDELTAG(ad,thetas(10))) + kDELTAG(ad,thetas(11));

modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 thetas([4 5])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(14)+thetas(15)*fp1fp2) + kGEOG(ad,thetas([3 6]))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([7 8])) + kDELTAG(ad,thetas(9))) + kDELTAG(ad,1).*(thetas(10).^2 + kMatG(dt1t2,thetas(11:13)));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

sqrt((thetTG(1)*100)^2 + thetTG(2)^2) 5 1e4 % global non-GIS amplitude
50 5 1e4 % GIS amplitude
thetTG(12) 5 1e4 % non-GIS regional amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

thetTG(13) 1 10 % geographic length scale

thetTG(10) 0 100 % linear regional amplitude
thetTG(11) .5 15 %  linear regional length scale
thetTG(7) 0 100 % linear local amplitude

100 0 1e4 % local offset

10 0 1e4 % local Matern
1.3 1 3e3
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
