function [TGdata2,TGdata,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(targcoords,addlsites,thetL0,winlength,thinlength,minlength,optimizemode,psmsldir,gslfile,maxdist,noiseMask)

% [TGdata2,TGdata,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(targcoords,
%             [addlsites],[thetL0],[winlength],[thinlength],[minlength],[optimizemode],
%             psmsldir,gslfile,[maxdist],[noiseMask])
%
% Find tide gauge sites to include, based on length and proximity criteria, then
% fit GP model to them in order to interpolate and take running averge.
%
% Specifically, from the site of all tide gauge sites with length of at least
% minlength(3), identifies sites that fit one of the following criteria:
%      1. are longer than minlength(1)
%      2. are within maxdist degrees of one of the sites listed in targcoords
%         and are longer than minlength(2)
%      3. is the nearest tide gauge to one of the sites listed in targcoords
%      4. has a PSMSL ID listed in addlsites
%
% Data are then smoothed and thinned by
% 
%      [TGdata2,thetL,TGmodellocal] = GPSmoothTideGauges(TGdata,winlength, ...
%                             optimizemode,thinlength,thetL0);
%
% EXAMPLE:
%
%  [TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords, ...
%                           [],[],[],[],[150 75 20],1.0,psmsldir,'none');
% 
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 21:04:21 EST 2014
%

defval('minlength',[150 75 20]);  % minimum length for global curves, most purposes and minimum
defval('thetL0',[]);
defval('winlength',11);
defval('thinlength',winlength-1);
defval('optimizemode',1.1); % set to 1.0 for local optimization only
defval('maxdist',5);
defval('thinyrstart',1700);
defval('noiseMask',[]);

defval('targcoords',[
34.97  -76.38; %Tump Point, NC
35.89  -75.64; %Sand Point, NC
39.085 -74.811; % Cape May Courthouse, NJ
39.50  -74.415; % Leeds Point, NJ
30.6   -81.7; % Nassau, FL
44.72  -63.23; % Chezzetcook, Nova Scotia
]);

defval('addlsites',[]);

defval('psmsldir',fullfile('.','rlr_annual'));
defval('gslfile',fullfile('.','CSIRO_Recons_gmsl_yr_2011.csv'));

% read PSMSL data for all sites longer than minlength(3)
[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(0,1000,minlength(3),psmsldir,gslfile);

% find long tide gauge sites anywhere
sub1=find((sitelen>minlength(1)));
sub1=union(sub1,find(TGsiteid==0));

% compute angular distances between sites, and add any tide gauges within 
% maxdist degrees with length greater than minlength(2)

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(mod(x1(:,2),360),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(mod(x2(:,2),360)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

sub2=[];
for ii=1:size(targcoords,1)
    TGdist=dDist(TGsitecoords,targcoords(ii,:));
    [m,mi]=min(TGdist);
    mi=union(mi,find((TGdist<maxdist).*(sitelen'>minlength(2))));
    sub2=union(sub2,mi);
end


% add addlsites
addlsites=addlsites(:)';
subaddlsites=find(ismember(TGsiteid,addlsites));

sitesub=union(sub1,sub2);
sitesub=union(sitesub,subaddlsites);
sub=find(ismember(TGid,TGsiteid(sitesub)));

clear TGdata;
TGdata.datid=TGid(sub);
TGdata.time1=TGcoords(sub,3);
TGdata.time2=TGcoords(sub,3);
TGdata.meantime=TGcoords(sub,3);
TGdata.limiting=zeros(size(TGid(sub)));
TGdata.Y=TGrsl(sub);
TGdata.dY =TGrslunc(sub);
TGdata.compactcorr=zeros(size(TGid(sub)));
TGdata.istg = ones(size(TGid(sub)));
TGdata.lat=TGcoords(sub,1);
TGdata.long=TGcoords(sub,2);;
TGdata.Ycv=sparse(diag(TGrslunc(sub).^2));
TGdata.siteid=TGsiteid(sitesub);
TGdata.sitenames=sitenames(sitesub);
TGdata.sitecoords=TGsitecoords(sitesub,:);
TGdata.sitelen=sitelen(sitesub);

%%%%%%%%%%

[TGdata2,thetL,TGmodellocal] = GPSmoothTideGauges(TGdata,winlength,optimizemode,thinlength,thetL0,thinyrstart,noiseMask);