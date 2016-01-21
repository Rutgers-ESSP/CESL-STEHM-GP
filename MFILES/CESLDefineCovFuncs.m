% This script defines some common covariance functions.
%
% kMat1 - Matern 1/2
% kMat3 - Matern 3/2
% kMat5 - Matern 5/2
% kSE - squared exponential
% kDP - dot product
% kCONST - constant
% kSIN - square exponential of sine function
% kMODSIN - product of kSIN and kSE
% kFREQ - cosine function
% kRQ - rational quadratic
% kMatG - Matern generic
% kDELTAG - delta function on a sphere
% kGEOG - Matern 5/2 on a sphere
% kGEOGG - Matern generic on a sphere
% 
% It also defines first and second derivatives kMat3d, kMat3dd and kDPd, kDPdd
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 28 16:09:12 EST 2015

% we will have a vector of observations
% and a vector of their uncertainty

% define covariance functions we will use

defval('refyear',1970);

dYears = @(years1,years2) abs(bsxfun(@minus,years1',years2));
dYears0 = @(years1,years2) (bsxfun(@minus,years1',years2));

angd = @(Lat0,Long0,lat,long) (180/pi)*(atan2( sqrt( (cosd(lat) .* sind(long-Long0)).^2 + (cosd(Lat0) .* sind(lat) - sind(Lat0) .* cosd(lat) .* cosd(long-Long0)).^2),(sind(Lat0) .* sind(lat) + cosd(Lat0) .* cosd(lat) .* cosd(long-Long0))));

dDist = @(x1,x2) angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))' + 1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

kMat1 = @(dx,thetas) thetas(1).^2 .* (1).*exp(-dx/thetas(2));

kMat3 = @(dx,thetas) thetas(1).^2 .* (1 + sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));

kMat5 = @(dx,thetas) thetas(1).^2 .* (1 + (sqrt(5)*dx/thetas(2)).*(1 + sqrt(5)*dx/thetas(2)/3)).*exp(-sqrt(5)*dx/thetas(2));

kSE = @(dx,thetas) thetas(1).^2 * exp(-(dx.^2)/(2*thetas(2).^2));

kDELTA = @(dx,thetas) thetas(1).^2 .* (dx==0);

kDP = @(years1,years2,thetas) thetas(1).^2 * bsxfun(@times,(years1-refyear)',(years2-refyear));
kCONST = @(thetas) thetas(1).^2;
kSIN = @(dt,thetas) thetas(1).^2 * exp(-2*sin(pi*dt/(thetas(2))).^2/thetas(3).^2);
kMODSIN = @(dt,thetas) kSIN(dt,thetas(1:3)) .* kSE(dt,[1 thetas(4)*thetas(2)]);
kFREQ = @(dt,thetas) thetas(1).^2 * cos(2*pi*dt/thetas(2));
kRQ = @(dx,thetas) thetas(1).^2 * (1 + dx.^2/(2*thetas(2)*thetas(3))).^-thetas(3);
kMatG = @(dx,thetas) thetas(1).^2 .* 2.^(1-thetas(3))./gamma(thetas(3)) .* (sqrt(2*thetas(3))*(dx+eps)/thetas(2)).^thetas(3) .* besselk(thetas(3),sqrt(2*thetas(3))*(dx+eps)/thetas(2));

kDELTAG = @(ad,thetas)thetas(1).^2.*(abs(ad)<1e-4).*(ad<360);

kGEOG = @(ad,thetas) kMat5(ad,thetas) .* (ad<360);
kGEOGG = @(ad,thetas)kMatG(ad,thetas).*(ad<360);


kMat3d = @(years1,years2,dx,thetas)thetas(1).^2.*(-3/(thetas(2).^2)).*dx.*exp(-sqrt(3)*dx/thetas(2)).*(-1+2*bsxfun(@ge,years1',years2));
kDPd = @(years1,years2,thetas) thetas(1).^2 * repmat((years1-refyear)',length(years2),1);
kMat3dd = @(dx,thetas)thetas(1).^2.*(3/(thetas(2).^2)).*(1-sqrt(3)*dx/thetas(2)).*exp(-sqrt(3)*dx/thetas(2));
kDPdd = @(years1,years2,thetas) thetas(1).^2 * ones(length(years2),length(years1));



FiniteMask = @(x1,x2) (bsxfun(@plus,sum(abs(x1),2)',sum(abs(x2),2)))<1e12;
