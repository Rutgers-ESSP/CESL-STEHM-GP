% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 17 19:22:54 EST 2014

% define some covariance functions

TGDefineCovFuncs;
defval('thetTG',[1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826]);


refyear=2010;
kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));

% asymptotic covariance associated with
% A exp((t0-t)/tau)
% so covariance of A^2 exp((2*t0-t1-t2)/tau)
%kASYMPT=@(years1,years2,thetas) thetas(1).^2*exp((2*thetas(2)-bsxfun(@plus,years1,years2'))/thetas(3));


modelspec(1).cvfunc = @(t1,t2,dt1t2,thetas,ad,fp1fp2) kMatG(dt1t2,[1 thetas([5 6])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(12)+thetas(13)*fp1fp2) + kGEOG(ad,thetas([3 7])) + kDELTAG(ad,thetas(4))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([8 9])) + kDELTAG(ad,thetas(10))) + kDELTAG(ad,thetas(11));

modelspec(1).traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) modelspec(1).cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluTGG = [

sqrt((thetTG(1)*100)^2 + thetTG(2)^2) 5 1e4 % global non-GIS amplitude
50 5 1e4 % GIS amplitude
thetTG(12) 5 1e4 % non-GIS regional amplitude
thetTG(8) 0 1e4 % local amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

thetTG(13) 1 10 % geographic length scale

%1e4 0 1e7 % asymptotic amplitude
%2500 2000 10000 % asymptotic end time
%7000 1000 50e3 % asymptotic timescale

thetTG(10) 0 100 % linear regional amplitude
thetTG(11) .5 15 % asympotic and linear regional length scale
thetTG(7) 0 100 % linear local amplitude

100 0 1e4 % local offset

1 1 1; % GIS global switch
1 1 1; % GIS regional switch

];

modelspec(1).thet0=tluTGG(:,1)';
modelspec(1).lb = tluTGG(:,2)';
modelspec(1).ub = tluTGG(:,3)';
modelspec(1).subfixed=[length(modelspec(1).thet0)-[1 0]]; % fix GIS switches
modelspec(1).sublength=[7 9];
