% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Dec 20 21:16:46 EST 2013

thinyrs=10;
minlen=50;

latNC1 =34.97; longNC1=-76.38; %Tump Point, NC
latNC2 = 35.89; longNC2=-75.64; %Sand Point, NC
latNJ1 = 39.085; longNJ1 = -74.811; % Cape May Courthouse, NJ
latNJ2 = 39.50; longNJ2 = -74.415; % Leeds Point, NJ

latNJ = mean([latNJ1 latNJ2]); longNJ = mean([longNJ1 longNJ2]);
latNC = mean([latNC1 latNC2]); longNC = mean([longNC1 longNC2]);
latFL = 30.6; longFL = -81.7;
latNS = 44.72; longNS = -63.23; % Chezzetcook, Nova Scotia


idNJ = 1e4;
idNC= 2e4;
idFL = 3e4;
idNS = 4e4;
idHolo = 2e5;

% First load tide gauge data

addlsites=[96 195 427 392 393]; % addl long Nova Scotia records
[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(960,960,[],[],[],addlsites);
TGsitelat=TGsitecoords(:,1);

notbedrock = union(TGsiteid(find(TGsitelat<39.266)),[1654 180 366 1637 519 875 848 362 856 430 351 776 367 1111 775]);
notbedrock = union(notbedrock,[idNJ idNC idNC+1 idNC+2 (idHolo+(7:16)*1e3) idFL ]);

notbedrock=setdiff(notbedrock,368);
bedrocksiteids=setdiff(TGsiteid,notbedrock);

bedrockMask = @(r1,r2) [repmat(sparse(~ismember(r1,bedrocksiteids)),1,length(r2)).*(repmat(sparse(~ismember(r2,bedrocksiteids))',length(r1),1))]';

% process tide gauge data

% 18 Apr: tide gauge
thetTG = [1.6585 36.1312 1e3 0.05 0 1.0217 0.5970   28.6848   11.9728    1.2888    5.7005   26.6414    4.1376   28.2188    7.3826];

%exclsites=[447 126 144 201 951 192];
exclsites=[];

[TGcoords,TGrsl,TGrslcv,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen,cvfuncTG] = ReadDenoisedPSMSLData(960,960,thetTG,bedrocksiteids,thinyrs,minlen,[],[],addlsites,exclsites);
TGyears=TGcoords(:,3);
TGlat=TGcoords(:,1);
TGlong=TGcoords(:,2);

% set up data structure


datid=[]; time1=[]; time2=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[];

datid = [datid ; TGid];
time1 = [time1 ; TGyears];
time2 = [time2 ; TGyears];
limiting = [limiting ; 0*TGrsl];
Y = [Y ; TGrsl];
dY = [dY ; TGrslunc];
compactcorr = [compactcorr ; 0*TGrsl];
istg = [istg ; ones(size(TGrsl))];
lat = [lat ; TGlat];
long = [long ; TGlong];

% process geological data

% NJ


datNJ=importdata(fullfile(IFILES,'Horton2013NJ.csv'));
dat2 = datNJ;

sub = find(dat2.data(:,2) ~= 40); % exclude Little Egg
wdatid = ones(length(sub),1)*idNJ;
count=[1:length(wdatid)]';
wtime1 = 1950-dat2.data(sub,8) + count/1e5;
wtime2 = 1950-dat2.data(sub,9) + count/1e5;
wlimiting = dat2.data(sub,3);
wY = dat2.data(sub,12)*1000;
wdY = dat2.data(sub,13)/2*1000;
wcompactcorr = dat2.data(sub,11)*1000;
wlat = dat2.data(sub,5);
wlong = dat2.data(sub,4);

datid = [datid ; wdatid];
time1 = [time1 ; wtime1];
time2 = [time2 ; wtime2];
limiting = [limiting ; wlimiting];
Y = [Y ; wY];
dY = [dY ; wdY];
compactcorr = [compactcorr ; wcompactcorr];
istg = [istg ; 0 * wY];
lat = [lat ; wlat];
long = [long ; wlong];

% NJ hi-res

datNJ2=importdata(fullfile(IFILES,'NJ_July2013.csv'));
dat2 = datNJ2;

wdatid=idNJ*ones(size(dat2.data(:,1)))+dat2.data(:,5);
count=[1:length(wdatid)]';
wtime1 = dat2.data(:,1) + dat2.data(:,2) + count/1e5;
wtime2 = dat2.data(:,1) - dat2.data(:,2) + count/1e5;
wY = dat2.data(:,3)*1000;
wdY = dat2.data(:,4)*1000;
wlimiting = 0*wY;
wcompactcorr = 0*wY*1000;
wlat = ones(size(wY)) * latNJ;
wlong = ones(size(wY)) * longNJ;

%  split out LP and CMCt

sub=find(wdatid == (idNJ+1));
wlat(sub) = latNJ1; wlong(sub)=longNJ1;
sub=find(wdatid == (idNJ+2));
wlat(sub) = latNJ2; wlong(sub)=longNJ2;

datid = [datid ; wdatid];
time1 = [time1 ; wtime1];
time2 = [time2 ; wtime2];
limiting = [limiting ; wlimiting];
Y = [Y ; wY];
dY = [dY ; wdY];
compactcorr = [compactcorr ; wcompactcorr];
istg = [istg ; 0 * wY];
lat = [lat ; wlat];
long = [long ; wlong];


% NC

datNC=importdata(fullfile(IFILES,'Kemp2011NCb.csv'));
dat2 = datNC;

wdatid=idNC*ones(size(dat2.data(:,1)))+dat2.data(:,1);
count=[1:length(wdatid)]';
wtime1 = dat2.data(:,2) + dat2.data(:,4) + count/1e5;
wtime2 = dat2.data(:,2) - dat2.data(:,4) + count/1e5;
wY = dat2.data(:,3)*1000;
wdY = dat2.data(:,5)*1000;
wlimiting = 0*wY;
wcompactcorr = 0*wY*1000;
wlat = ones(size(wY)) * latNC;
wlong = ones(size(wY)) * longNC;

%  split out Tump Point and Sand Point

sub=find(wdatid == idNC+1);
wlat(sub) = latNC1; wlong(sub)=longNC1;
sub=find(wdatid == idNC+2);
wlat(sub) = latNC2; wlong(sub)=longNC2;

datid = [datid ; wdatid];
time1 = [time1 ; wtime1];
time2 = [time2 ; wtime2];
limiting = [limiting ; wlimiting];
Y = [Y ; wY];
dY = [dY ; wdY];
compactcorr = [compactcorr ; wcompactcorr];
istg = [istg ; 0 * wY];
lat = [lat ; wlat];
long = [long ; wlong];

% FL

datFL=importdata(fullfile(IFILES,'NassauFL_July2013.csv'));
dat2 = datFL;

wdatid=idFL*ones(size(dat2.data(:,1)));
count=[1:length(wdatid)]';
wtime1 = dat2.data(:,1) + dat2.data(:,2) + count/1e5;
wtime2 = dat2.data(:,1) - dat2.data(:,2) + count/1e5;
wY = dat2.data(:,3)*1000;
wdY = dat2.data(:,4)*1000;
wlimiting = 0*wY;
wcompactcorr = 0*wY*1000;
wlat = ones(size(wY)) * latFL;
wlong = ones(size(wY)) * longFL;

datid = [datid ; wdatid];
time1 = [time1 ; wtime1];
time2 = [time2 ; wtime2];
limiting = [limiting ; wlimiting];
Y = [Y ; wY];
dY = [dY ; wdY];
compactcorr = [compactcorr ; wcompactcorr];
istg = [istg ; 0 * wY];
lat = [lat ; wlat];
long = [long ; wlong];

% NS Gehrels et al. 2005

datNS=importdata(fullfile(IFILES,'Gehrels2005_NS.csv'));
dat2 = datNS;

wdatid=idNS*ones(size(dat2.data(:,1)));
count=[1:length(wdatid)]';
wtime1 = dat2.data(:,1) + dat2.data(:,2) + count/1e5;
wtime2 = dat2.data(:,1) - dat2.data(:,3) + count/1e5;
wY = dat2.data(:,4)*1000;
wdY = dat2.data(:,5)*1000;
wlimiting = 0*wY;
wcompactcorr = 0*wY*1000;
wlat = ones(size(wY)) * latNS;
wlong = ones(size(wY)) * longNS;

datid = [datid ; wdatid];
time1 = [time1 ; wtime1];
time2 = [time2 ; wtime2];
limiting = [limiting ; wlimiting];
Y = [Y ; wY];
dY = [dY ; wdY];
compactcorr = [compactcorr ; wcompactcorr];
istg = [istg ; 0 * wY];
lat = [lat ; wlat];
long = [long ; wlong];

% Engelhart & Horton database

datHolo=importdata(fullfile(IFILES,'Engelhart_Horton_2012.csv'));
HoloRegions=[1:7 9:16];
for curreg=1:length(HoloRegions)
	sub=find((datHolo.data(:,1)==HoloRegions(curreg)).*(datHolo.data(:,2)==0));

	count=[1:length(sub)]';
	wdatid = ones(length(sub),1)*(idHolo+HoloRegions(curreg)*1e3);
	wtime1=1950-datHolo.data(sub,8) + count/1e5;
	wtime2=1950-datHolo.data(sub,9) + count/1e5;
	wlimiting=datHolo.data(sub,2);
	wY=datHolo.data(sub,10)*1000;
	wdY=datHolo.data(sub,11)*1000;
	wlat=datHolo.data(sub,3);
	wlong=-datHolo.data(sub,4);
	wcompactcorr = wY*.1;
	
	datid = [datid ; wdatid];
	time1 = [time1 ; wtime1];
	time2 = [time2 ; wtime2];
	limiting = [limiting ; wlimiting];
	Y = [Y ; wY];
	dY = [dY ; wdY];
	compactcorr = [compactcorr ; wcompactcorr];
	istg = [istg ; 0 * wY];
	lat = [lat ; wlat];
	long = [long ; wlong];
end

meantime=mean([time1 time2],2);
dt = abs(time1-time2)/4;

Ycv = sparse(diag(dY.^2));
Ycv(1:length(TGrsl),1:length(TGrsl)) = TGrslcv;

dY0=dY;
Ycv0 = Ycv;

compactcorr=sparse(compactcorr);

% subtract GIA model

giamodel.gia=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Dsea_250');
giamodel.lat=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lat');
giamodel.long=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lon');

ICE5Gin=giamodel.gia;
ICE5Glat=giamodel.lat;
ICE5Glon=giamodel.long;
sub=find(ICE5Glon>180); ICE5Glon(sub)=ICE5Glon(sub)-360;
[ICE5Glon,si]=sort(ICE5Glon); ICE5Gin=ICE5Gin(si,:);

[regionsu,regionsusi]=unique(datid);
sitecoords=[lat(regionsusi) long(regionsusi)];
GIAproju=zeros(size(regionsu));
GIAproj=zeros(size(Y));
for i=1:length(GIAproju)
	if regionsu(i)>0
		GIAproju(i)=interp2(ICE5Glat,ICE5Glon,ICE5Gin,sitecoords(i,1),sitecoords(i,2));
		sub=find(datid==regionsu(i));
		GIAproj(sub)=GIAproju(i).*(meantime(sub)-1970);
	end
end
Y0=Y;
Y=Y0-GIAproj;

% indices

datNAO = importdata(fullfile(IFILES,'nao-trouet2009.txt'));
	
	%% comparison to temperature records
	impt=importdata(fullfile(IFILES,'Marcott2013_global.txt'));
	Marcottgl_yr = 1950-impt.data(:,1);
	Marcottgl_T = impt.data(:,2);
	Marcottgl_dT = impt.data(:,3);

		impt=importdata(fullfile(IFILES,'Marcott2013_regional.txt'));
	Marcottr_yr=1950-impt.data(:,1);
	Marcott_Next_T=impt.data(:,2);
	Marcott_Next_dT=impt.data(:,3);
	Marcott_eq_T=impt.data(:,4);
	Marcott_eq_dT=impt.data(:,5);
	Marcott_Sext_T=impt.data(:,6);
	Marcott_Sext_dT=impt.data(:,7);

	impt=importdata(fullfile(IFILES,'HadCRUT4-gl.dat'));
	HadCRUT_yr=impt(1:2:end,1);
	HadCRUT_T=impt(1:2:end,end);

	Mann_T=importdata(fullfile(IFILES,'glglfulihad_smxx.txt'));
	Mann_yr = [[1:length(Mann_T)]-1]';

% fingerprint

[GISfp,GISfplong,GISfplat]=readFingerprintInd('gis',IFILES);
sub=find(GISfplong>180); GISfplong(sub)=GISfplong(sub)-360;
[GISfplong,si]=sort(GISfplong); GISfp=GISfp(:,si);
GISfp=GISfp*1000;


obsGISfp = interp2(GISfplong,GISfplat,GISfp,long,lat,'linear');
obsGISfp(find(lat>100))=1;
