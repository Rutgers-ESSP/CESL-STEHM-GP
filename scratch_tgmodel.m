% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 24 11:25:42 EST 2014
%
% Find tide gauge sites to include, based on length and proximity criteria, then
% fit GP model to them in order to interpolate and take running averge.
%
% Should tweak to use proxy data set rather than prescribed list of sites to look for 
% tide gauges near.

ROOTDIR='~/Dropbox/Consulting/Risky Business/Code/RiskyBusinessScience/slr';
addpath(fullfile(ROOTDIR,'MFILES'));
IFILES=fullfile(ROOTDIR,'../../IFILES/slr/');
PARAMDIR=fullfile(ROOTDIR,'PARAMS/');
OUTDIR=[pwd '/'];
GIAanchoryear=2005;

giafile=fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc');
psmsldir=fullfile(IFILES,'rlr_annual');
gslfile=fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv');

giamodel.gia=ncread(giafile,'Dsea_250');
giamodel.lat=ncread(giafile,'Lat');
giamodel.long=ncread(giafile,'Lon');

minlength=70; % minimum length for most purposes
minlength0=20; % absolute minimum length
thetG=[1.621 11.966 40.225 0.549];

[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(1,1000,minlength0,psmsldir,gslfile);
sub1=find((sitelen>130));
sub1=union(sub1,find(TGsiteid==0));

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

targcoords=[
34.97  -76.38; %Tump Point, NC
35.89  -75.64; %Sand Point, NC
39.085 -74.811; % Cape May Courthouse, NJ
39.50  -74.415; % Leeds Point, NJ
30.6   -81.7; % Nassau, FL
44.72  -63.23; % Chezzetcook, Nova Scotia
];

sub2=[];
for ii=1:size(targcoords,1)
    TGdist=dDist(TGsitecoords,targcoords(ii,:));
    [m,mi]=min(TGdist);
    mi=union(mi,find((TGdist<2).*(sitelen'>minlength)));
    sub2=union(sub2,mi);
end

addlsites=[
350 % Portsmouth, UK
429 % New London, CT
1453 % Rarotonga, Cook Islands
638 % Reyjkavik
503 % Alexandria, Egypt
440 % Eugene Island, LA
183 % Portland, ME
235 % Boston, MA
150 % Auckland, NZ
826 % Simons Bay, South Africa
488 % Tarifa, Spain
65  % Sydney, Australia
1216 % Spring Bay, Tasmania
983 % Cocos Island
];
addlsites=addlsites';
subaddlsites=find(ismember(TGsiteid,addlsites));

sitesub=union(sub1,sub2);
sitesub=union(sitesub,subaddlsites);
sub=find(ismember(TGid,TGsiteid(sitesub)));

clear TGdata;
TGdata.datid=TGid(sub);
TGdata.time1=TGcoords(sub,3);
TGdata.time2=TGcoords(sub,3);
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

dospatial=0;
Nlongest=[20];

[thetGLR,cvfuncGLR,thetGLRA] = TrainGPSLModel([TGdata.lat TGdata.long TGdata.time1],TGdata.Y,TGdata.dY,TGdata.datid,TGdata.siteid,TGdata.sitenames,TGdata.sitecoords,TGdata.sitelen,minlength,giamodel,thetG,dospatial,Nlongest);

% 24 Feb 2014
%thetGLR=[1.621 11.966 40.225 0.549 1.540 7.607 0.060 38.927 0.938 0.500 0.264 6.548];

%%

clear TGmodel TGtestsitedef;

TGmodel.cvfunc=@(t1,t2,dt1t2,thetas,dy1y2,bedMsk,fp1fp2) cvfuncGLR(t1,t2,thetas,dy1y2,dt1t2);
TGmodel.traincv = @(t1,t2,dt1t2,thetas,errcv,ad,bedmask,fp1fp2) TGmodel.cvfunc(t1,t2,dt1t2,thetas,ad,bedmask,fp1fp2) + errcv;
TGmodel.thet0=thetGLR;

TGtestsitedef.sites=[TGdata.siteid TGdata.sitecoords];
TGtestsitedef.names=TGdata.sitenames;
TGtestsitedef.names2=TGtestsitedef.names;

winlength=11;
thinlength=10;
for ii=1:length(TGdata.siteid)
    sub=find(TGdata.datid==TGdata.siteid(ii));
    TGtestsitedef.t{ii} = floor(min(TGdata.time1(sub)-winlength)):ceil(max(TGdata.time1(sub)+winlength));
end
%% now do a regression

noiseMasks = ones(1,length(thetGLR));
noiseMasklabels={'full'};

matlabpool close;

trainsub=find(TGdata.limiting==0);
[TGf,TGsd,TGV,TGtestlocs]=RegressHoloceneDataSets(TGdata,TGtestsitedef,TGmodel,thetGLR,[],giamodel,trainsub,[],noiseMasks,TGtestsitedef.t,GIAanchoryear);

% now take running averages
Mop=abs(bsxfun(@minus,TGtestlocs.X(:,3),TGtestlocs.X(:,3)'))<(winlength/2);
Mop=Mop.*bsxfun(@eq,TGtestlocs.reg,TGtestlocs.reg');
adder=sum(Mop,2);
sub=find(adder==winlength);
Mop=Mop(sub,:)/winlength;
t2=Mop*TGtestlocs.X(:,3);
reg2=round(Mop*TGtestlocs.reg*10)/10;
selsub=[];
for ii=1:length(TGdata.siteid)
    subO=find(TGdata.datid==TGdata.siteid(ii));
    sub=find((reg2==TGdata.siteid(ii)).*(t2>=min(TGdata.time1(subO))).*(t2<=max(TGdata.time1(subO))));
    sub3=sub(1:thinlength:end);
    selsub=[selsub ; sub3];
end
Mop=Mop(selsub,:);

TGf2=Mop*TGf;
TGV2=Mop*TGV*Mop';
TGsd2=sqrt(diag(TGV2));

TGdata2.datid=Mop*TGtestlocs.reg;
TGdata2.time1=Mop*TGtestlocs.X(:,3);
TGdata2.time2=TGdata2.time1;
TGdata2.limiting=zeros(size(TGdata2.datid));
TGdata2.Y=TGf2;
TGdata2.dY =TGsd2;
TGdata2.compactcorr=zeros(size(TGdata2.datid));
TGdata2.istg = ones(size(TGdata2.datid));
TGdata2.lat=Mop*TGtestlocs.X(:,1);
TGdata2.long=Mop*TGtestlocs.X(:,2);
TGdata2.Ycv=sparse(TGV2);
TGdata2.siteid=TGdata.siteid;
TGdata2.sitenames=TGdata.sitenames;
TGdata2.sitecoords=TGdata.sitecoords;
TGdata2.sitelen=TGdata.sitelen;


    

