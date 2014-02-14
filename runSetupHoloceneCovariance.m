% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 11 18:56:57 EST 2014

% define some covariance functions

TGDefineCovFuncs;
bedmsk = bedrockMask(datid,datid);

% asymptotic covariance associated with
% A exp((t0-t)/tau)
% so covariance of A^2 exp((2*t0-t1-t2)/tau)
kASYMPT=@(years1,years2,thetas) thetas(1).^2*exp((2*thetas(2)-bsxfun(@plus,years1,years2'))/thetas(3));


%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask) cvfuncTG(t1,t2,dt1t2,[thetas(1:4) 0 1 thetas(5:6) 0 thetas(7:10) 0 1] ,ad,bedmask) + kMat3(dt1t2,thetas(11:12)) + kGEOG(ad,[1 thetas(15)]) .* kMat5(dt1t2,thetas(13:14)) + kGEOG(ad,[1 thetas(10)]).* bedmask .* kMat5(dt1t2,thetas(16:17));

%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask) cvfuncTG(t1,t2,dt1t2,[thetas(1:4) 0 1 thetas(5:6) 0 thetas(7:10) 0 1] ,ad,bedmask) + kMat3(dt1t2,thetas(11:12)) + kGEOG(ad,[1 thetas(16)]) .* kMatG(dt1t2,thetas(13:15));

%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask) kRQ(dt1t2,[1,thetas(2:3)]).*(thetas(1).^2+kDELTAG(ad,thetas(4))+kGEOG(ad,thetas(5:6))) +  kDP(t1,t2,1).*(kDELTAG(ad,thetas(7)).*bedmask+kGEOG(ad,thetas(8:9))) + kMat3(dt1t2,thetas(10:11)) + kGEOG(ad,[1 thetas(15)]) .* kRQ(dt1t2,thetas(12:14));

%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask) kMatG(dt1t2,[1 thetas(1:2)]).*(thetas(3).^2 + kGEOG(ad,[thetas(4:5)])) + kMatG(dt1t2,[1 thetas(6:7)]).*(thetas(8).^2 + kGEOG(ad,thetas([9 5])) + kDELTAG(ad,thetas(10))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([11 12])) + kDELTAG(ad,thetas(13)));


%tluTGG = [
%thetTG(1)   0    100 % DP
%
%thetTG(2)  0   1e4 %RQ
%thetTG(3)   5   1e3
%thetTG(4)  .05  10
%
%thetTG(7)    0   100   % local amplitudes
%thetTG(8)   0   1e4
%
%thetTG(10)    0   100	%DP geog amplitude and length
%thetTG(11)    .5 100
%
%thetTG(12)   0   1e4	%RQ geog amplitude and length
%thetTG(13)    1 100
%
%4e3 0   100e3 % Mat3 global
%3e3 100 30e3
%
%4e3 0 100e3 % MatG regional
%3e3 100 30e3
%2 .5 7.5
%thetTG(11)    1 100
%
%%1e2  0 20e3 % Mat5 regional (dynamics)
%%1e3 10 30e3
%
%];


%tluTGG = [
%thetTG(2)  0   1e4 %RQ
%thetTG(3)   5   1e3
%thetTG(4)  .05  10
%
%thetTG(8)   0   1e4 % RQ local
%
%thetTG(12)   0   1e4	%RQ geog amplitude and length
%thetTG(13)    1 100
%
%thetTG(7)    0   100   % DP local
%
%thetTG(10)    0   100	%DP geog amplitude and length
%thetTG(11)    .5 15
%
%4e3 0   100e3 % Mat3 global
%1e3 100 3e3
%
%4e3 10 100e3 % RQ regional
%3e3 100 30e3
%1   .01  10
%thetTG(11)    1 15
%
%];

%
%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask) thetas(1).^2.*(kMatG(dt1t2,thetas([2 3 4])) + kMatG(dt1t2,[(1-thetas(2)) thetas([5 6])])) + kGEOG(ad,thetas(7:8)).*(kMatG(dt1t2,thetas([9 [3 4]])) + kMatG(dt1t2,[(1-thetas(9)) thetas([5 6])])) + kDELTAG(ad,thetas(10)).*(kMatG(dt1t2,[1 thetas([5 6])])) + kDP(t1,t2,1).*(kGEOG(ad,thetas([11 12])) + bedmask.*kDELTAG(ad,thetas(13)));
%
%traincvTGG = @(t1,t2,dt1t2,thetas,errcv,ad,bedmask) cvfuncTGG(t1,t2,dt1t2,thetas,ad,bedmask) + errcv;
%
%tluTGG = [
%
%sqrt((thetTG(1)*100)^2 + thetTG(2)^2) 0 1e4 % global amplitude
%
%0.5 0 1 % global centennial fraction
%
%1e3 100 3e3 % centennial temporal parameters
%1.5 .5 9.5
%
%60 1 100 % decadal temporal parameters
%1.5 .5 9.5
%
%thetTG(12) 0 1e4 % regional amplitude
%thetTG(13) 1 100 % Matern regional length scale
%
%0.5 0 1 % regional centennial fraction
%
%thetTG(8) 0 1e4 % local decadal amplitude
%
%thetTG(10) 0 100 % linear regional amplitude
%thetTG(11) .5 15 % linear regional length scale
%thetTG(7) 0 100 % linear local amplitude
%]
%
%
%
%
%
%thetTGG=tluTGG(:,1)';
%lbTGG = tluTGG(:,2)';
%ubTGG = tluTGG(:,3)';
%%subfixedTGG = [1:10];
%%subfixedTGG=setdiff(subfixedTGG,[1 5 7]);
%%subfixedTGG=1:9;
%subfixedTGG=[8 10:13]; % fix GIA factors, local amplitude, regional length scale 

cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask,fp1fp2) kMatG(dt1t2,[1 thetas([5 6])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(12)+thetas(13)*fp1fp2) + kGEOG(ad,thetas([3 7])) + kDELTAG(ad,thetas(4))) + kDP(t1,t2,1).*(kGEOG(ad,thetas([8 10])) + bedmask.*kDELTAG(ad,thetas(9))) + kDELTAG(ad,thetas(11));

% TBD!!
%cvfuncTGG = @(t1,t2,dt1t2,thetas,ad,bedmask,fp1fp2) kMatG(dt1t2,[1 thetas([5 6])]) .*(thetas(1).^2 + thetas(2).^2*(thetas(12)+thetas(13)*fp1fp2) + kGEOG(ad,thetas([3 7])) + kDELTAG(ad,thetas(4))) +
%kASYMPT(t1,t2,thetas([8:10])).*kGEOG(ad,[1 thetas(11))
% kDP(t1,t2,1).*(kGEOG(ad,thetas([8 10])) + bedmask.*kDELTAG(ad,thetas(9))) + kDELTAG(ad,thetas(11));


traincvTGG = @(t1,t2,dt1t2,thetas,errcv,ad,bedmask,fp1fp2) cvfuncTGG(t1,t2,dt1t2,thetas,ad,bedmask,fp1fp2) + errcv;

tluTGG = [

sqrt((thetTG(1)*100)^2 + thetTG(2)^2) 5 1e4 % global non-GIS amplitude
50 5 1e4 % GIS amplitude
thetTG(12) 5 1e4 % non-GIS regional amplitude
thetTG(8) 0 1e4 % local amplitude

1e3 1 3e3 % temporal parameters
1.5 .5 9.5

thetTG(13) 1 10 % geographic length scale

thetTG(10) 0 100 % linear regional amplitude
thetTG(7) 0 100 % linear local amplitude
thetTG(11) .5 15 % linear regional length scale

100 0 1e4 % local offset

1 1 1; % GIS global switch
1 1 1; % GIS regional switch


]





thetTGG=tluTGG(:,1)';
lbTGG = tluTGG(:,2)';
ubTGG = tluTGG(:,3)';
%subfixedTGG = [1:10];
%subfixedTGG=setdiff(subfixedTGG,[1 5 7]);
%subfixedTGG=1:9;
%subfixedTGG=[ 8:9]; % fix GIA length scale 
subfixedTGG=[12 13]; % fix GIS switches
