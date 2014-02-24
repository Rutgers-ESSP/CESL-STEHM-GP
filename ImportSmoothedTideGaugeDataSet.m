function [TG,cvfuncTG,notbedrock]=ImportSmoothedTideGaugeDataSet(includeCW,thinyrs,minlen)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 17 23:16:09 EST 2014


defval('includeCW',0);
defval('thinyrs',10);
defval('minlen',50);

% set up bedrockMask
addlsites=[96 195 427 392 393]; % addl long Nova Scotia records
[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(960,960,[],[],[],addlsites);
TGsitelat=TGsitecoords(:,1);

notbedrock = union(TGsiteid(find(TGsitelat<39.266)),[1654 180 366 1637 519 875 848 362 856 430 351 776 367 1111 775]);
notbedrock=setdiff(notbedrock,368);
bedrocksiteids=setdiff(TGsiteid,notbedrock);

% process tide gauge data

% 18 Apr: tide gauge
thetTG = [1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826];

%exclsites=[447 126 144 201 951 192];
exclsites=[];
if ~includeCW
    exclsites=[exclsites 0];
end

[TGcoords,TG.Y,TG.Ycv,TG.dY,TG.datid,TG.siteid,TG.sitenames,TG.sitecoords,TG.sitelen,cvfuncTG] = ReadDenoisedPSMSLData(960,960,thetTG,bedrocksiteids,thinyrs,minlen,[],[],addlsites,exclsites);
TG.years=TGcoords(:,3);
TG.lat=TGcoords(:,1);
TG.long=TGcoords(:,2);