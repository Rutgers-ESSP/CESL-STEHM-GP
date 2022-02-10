% Read in data files
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:06:25 -0400

minperstudy=1; % minimum size of studies to include

loadNearbyTG=1; % load nearby tide gauge sites?
psmsldir=fullfile(IFILES2,'rlr_annual'); % location of PSMSL data files
thinyrs=10; % how much to thin/smooth tide gauge data
optimizemode=1.0; % local optimization only for smoothing GPs for tide gauges

usePriorGIA=1; % use a GIA model as the mean of the prior?
refyear=2010; % reference year for GIA hindcast

datid=[]; time1=[]; time2=[]; mediantime=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[];
siteid=[]; sitenames={}; sitecoords=[];

% import data file
% this may require some adjustment for different applications

datPX = importdata(PXdatafile);
datPX.textdata=datPX.textdata(2:end,:);

% catch entries without age errors
sub=find(isnan(datPX.data(:,8))); datPX.data(sub,8)=10000;
sub=find(isnan(datPX.data(:,7))); datPX.data(sub,7)=10000;

study=datPX.textdata(:,1);
for www=1:length(study)
    txtsub=strfind(study{www},',');
    if length(txtsub)>0
        study{www}=study{www}(1:txtsub(1)-1);
    end
end

uStudy = unique(study);
for ii=1:length(uStudy)
    sub=find(strcmpi(uStudy{ii},study));
    if length(sub)>=minperstudy
        site = datPX.textdata(sub,2);
        uSite=unique(site);
        for jj=1:length(uSite)
            curid = 1e4*ii + jj;
            curstudysite=[uStudy{ii} '-' uSite{jj}];
            sub2=sub(find(strcmpi(uSite{jj},site)));
            wdatid=ones(length(sub2),1)*curid;
            wmediantime=datPX.data(sub2,6);
            wtime1=wmediantime-datPX.data(sub2,8)+(1:length(sub2))'/1e5;
            wtime2=wmediantime+datPX.data(sub2,7)+(1:length(sub2))'/1e5;
            
            wlimiting=zeros(length(sub2),1);
            wYmedian=datPX.data(sub2,3);
            wY1=wYmedian-datPX.data(sub2,5);
            wY2=wYmedian+datPX.data(sub2,4);
            wY=(wY1+wY2)/2;
            wdY=abs(wY2-wY1)/2;
            wcompactcorr=zeros(length(sub2),1);;
            wlat=datPX.data(sub2,1);
            wlong=datPX.data(sub2,2);
            
            wY=wY*1000;
            wdY=wdY*1000;
            
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
            mediantime = [mediantime ; wmediantime];

            siteid=[siteid ; curid];
            sitenames={sitenames{:}, curstudysite};
            sitecoords=[sitecoords; mean(wlat) mean(wlong)];
        end

    end
end 



PX.datid=round(datid);
PX.time1=time1;
PX.time2=time2;
PX.limiting=limiting;
PX.Y=Y;
PX.dY = dY;
PX.compactcorr=compactcorr;
PX.istg = istg;
PX.lat=lat;
PX.long=long;
PX.Ycv=sparse(diag(dY.^2));
PX.siteid=round(siteid);
PX.sitenames=sitenames;
PX.meantime=(PX.time1+PX.time2)/2;
PX.sitecoords = sitecoords;

PX0=PX;

% if you want to subset in any way, do it here
% currently just subsetting those with latitude <= 90 (i.e., everything)
sub=find((PX.lat<=latlim(2)).*(PX.lat>=latlim(1)).*(PX.long>=longlim(1)).*(PX.long<=longlim(2)));
subS=find((PX.sitecoords(:,1)<=latlim(2)).*(PX.sitecoords(:,1)>=latlim(1)).*(PX.sitecoords(:,2)>=longlim(1)).*(PX.sitecoords(:,2)<=longlim(2)));

PX=SubsetDataStructure(PX,sub,subS);


%old tide gauges

clear TGold;
TGold.datid = [];
TGold.meantime=[];
TGold.Y=[];
TGold.dY=[];
TGold.lat=[];
TGold.long=[];
TGold.siteid=[];
TGold.sitenames={};
TGold.sitecoords=[];
TGold.sitelen=[];

for ppp=1:3
    if ppp==1
        dat=importdata(fullfile(IFILES2,'amsterdam.sea.level.txt'));
        wsite='AMSTERDAM_OLD'; wcurid=5001;
        wlat = [52.3667]; wlong=[4.9000];
        wdatid = ones(length(dat.data),1)*wcurid;
        wtime = dat.data(:,1);
        wY = dat.data(:,2);
    elseif ppp==2
        dat=importdata(fullfile(IFILES2,'Kronstadt_ReportsFGI_Bogdanov_appendix.csv'));
        wsite='KRONSTADT_OLD'; wcurid=5002;
        wlat = [59.98]; wlong=[29.77];
        wdatid = ones(length(dat.data),1)*wcurid;
        wtime = dat.data(:,1);
        wY = dat.data(:,3);
    elseif ppp==3
         dat=importdata(fullfile(IFILES2,'ekman_2003_stockholm.csv'));
        wsite='STOCKHOLM_OLD'; wcurid=5003;
        wlat = [59.32]; wlong=[18.08];
        wdatid = ones(length(dat.data),1)*wcurid;
        wtime = dat.data(:,1);
        wY = dat.data(:,3);
    end
    
    TGold.datid = [TGold.datid;wdatid];
    TGold.meantime = [TGold.meantime; wtime];
    TGold.Y=[TGold.Y ; wY];
    TGold.dY = ones(size(TGold.Y))*3;
    TGold.lat = [TGold.lat ; ones(size(wY))*wlat];
    TGold.long = [TGold.lat ; ones(size(wY))*wlong];
    TGold.siteid = [TGold.siteid ; wcurid];
    TGold.sitenames={TGold.sitenames{:},wsite};
    TGold.sitecoords=[TGold.sitecoords ; wlat wlong];
    TGold.sitelen = [TGold.sitelen ; length(wY)]; 

 end

TGold.time1 = TGold.meantime;
TGold.time2 = TGold.meantime;
TGold.limiting = zeros(size(TGold.datid));
TGold.compactcorr = zeros(size(TGold.datid));
TGold.istg = ones(size(TGold.datid));
TGold.Ycv=sparse(diag(TGold.dY.^2));

TGnoisemask = [];
TGoldS = GPSmoothTideGauges(TGold,11,1.0,10,[],1700,TGnoisemask);


TGnoisemask = [];
if loadNearbyTG
    [TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords,[],[],[],[],[150 75 20],optimizemode,psmsldir,'none',[],TGnoisemask);
end
    
% drop Trois-Rivieres
sub=find(TG.datid~=126);
subS=find(TG.siteid~=126);
TG=SubsetDataStructure(TG,sub,subS);
TG0=SubsetDataStructure(TG0,sub,subS);


% Hay et al 2015 GMSL curve

gslfile=fullfile(IFILES2,'Hay2015_KFandGP_GMSL.mat');

Haydat=load(gslfile);
Hay.Y=Haydat.KF_GMSL(:);
Hay.dY=sqrt(diag(Haydat.KF_GMSL_var));
Hay.Ycv=Haydat.KF_GMSL_var;
Hay.datid=0*Hay.Y;
Hay.time1=Haydat.tt_KF;
Hay.time2=Hay.time1;
Hay.meantime=Hay.time1;
Hay.lat=ones(size(Hay.Y))*1e6;
Hay.long=Hay.lat;
Hay.compactcorr=0*Hay.Y;
Hay.limiting=0*Hay.Y;
Hay.istg=ones(size(Hay.Y));
Hay.siteid=0;
Hay.sitenames={'Hay_KF_GMSL'};
Hay.sitecoords=[1e6 1e6];
Hay.sitelen=length(Hay.Y);

% smooth Hay et al curve
Hayavgwin=10;
Haystep=10;
HayGSL=Hay;
HayGSL.time1=[1885:Haystep:2005]';
HayGSL.time2=HayGSL.time1; HayGSL.meantime=HayGSL.time1;
M=abs(bsxfun(@minus,HayGSL.time1,Hay.time1))<=(Hayavgwin/2);
M=bsxfun(@rdivide,M,sum(M,2));
HayGSL.Y=M*Hay.Y;
HayGSL.Ycv=M*Hay.Ycv*M';
HayGSL.dY=sqrt(diag(HayGSL.Ycv));
HayGSL.datid=0*HayGSL.Y;
HayGSL.lat=ones(size(HayGSL.Y))*1e6;
HayGSL.long=HayGSL.lat;
HayGSL.compactcorr=0*HayGSL.Y;
HayGSL.limiting=0*HayGSL.Y;
HayGSL.istg=ones(size(HayGSL.Y));

% add in GSL flattener
% so that average rate of change between -100 to 100 CE and 1600-1800 CE is close to zero

clear GSLflattener;
GSLflattener.sigma=1e4;
GSLflattener.time1=.01+[-100:50:100 1600:50:1800]';
GSLflattener.Y=0*GSLflattener.time1;
GSLflattener.Ycv = ones(length(GSLflattener.Y))*(GSLflattener.sigma^2);
GSLflattener.Ycv(1:5,1:5) = GSLflattener.Ycv(1:5,1:5)+(.0025*1700)^2+eye(5)*GSLflattener.sigma^2-GSLflattener.sigma^2/5;
GSLflattener.Ycv(6:10,6:10) = GSLflattener.Ycv(6:10,6:10)+eye(5)*GSLflattener.sigma^2-GSLflattener.sigma^2/5;

GSLflattener.dY=sqrt(diag(GSLflattener.Ycv));
GSLflattener.datid=0*GSLflattener.Y;
GSLflattener.time2=GSLflattener.time1;
GSLflattener.meantime=GSLflattener.time1;
GSLflattener.lat=ones(size(GSLflattener.Y))*1e6;
GSLflattener.long=GSLflattener.lat;
GSLflattener.compactcorr=0*GSLflattener.Y;
GSLflattener.limiting=0*GSLflattener.Y;
GSLflattener.istg=ones(size(GSLflattener.Y));
GSLflattener.siteid=0;
GSLflattener.sitenames={'GSLflattener'};
GSLflattener.sitecoords=[1e6 1e6];
GSLflattener.sitelen=length(GSLflattener.Y);

% add in old tide gauges
TG=MergeDataStructures(TG,TGoldS);

% add in GMSL curve
TGNOGSL=TG;
TG=MergeDataStructures(TG,HayGSL);

%% load GIA model

if usePriorGIA
    giamodel.gia=ncread(fullfile(IFILES2,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Dsea_250');
    giamodel.lat=ncread(fullfile(IFILES2,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lat');
    giamodel.long=ncread(fullfile(IFILES2,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lon');

    ICE5Ggia=giamodel.gia;
    ICE5Glat=giamodel.lat;
    ICE5Glon=giamodel.long;
    sub=find(ICE5Glon>180); ICE5Glon(sub)=ICE5Glon(sub)-360;
    [ICE5Glon,si]=sort(ICE5Glon); ICE5Ggia=ICE5Ggia(si,:);
else
    ICE5Ggia=0;
    ICE5Glat=0;
    ICE5Glon=0;
end


% create data structures

clear datasets;

if loadNearbyTG
    datasets{1}=MergeDataStructures(TG,PX);
    datasets{2}=MergeDataStructures(MergeDataStructures(TG,PX),GSLflattener);
    datasets{3}=PX;
    datasets{4}=MergeDataStructures(TGNOGSL,PX);

    datasets{1}.label='TG+GSL+PX';
    datasets{2}.label='TG+GSL+PX+flat';
    datasets{3}.label='PX';
    datasets{4}.label='TG+PX';
else
    datasets{1}=PX;
    datasets{1}.label='PX';
end

for ii=1:length(datasets)
    t1=datasets{ii}.time1; t2=datasets{ii}.time2;
    datasets{ii}.long = mod(datasets{ii}.long,360); sub=find(datasets{ii}.long>180); datasets{ii}.long(sub)=datasets{ii}.long(sub)-360;
    datasets{ii}.meantime=mean([t1 t2],2);
    datasets{ii}.dt = abs(t1-t2)/4;
    datasets{ii}.compactcorr=sparse(datasets{ii}.compactcorr);

    if usePriorGIA
        % subtract GIA model
        datasets{ii}.siteGIA = interp2(ICE5Glat,ICE5Glon,ICE5Ggia, ...
                                       datasets{ii}.sitecoords(:,1),datasets{ii}.sitecoords(:,2));
        datasets{ii}.siteGIA(find(datasets{ii}.sitecoords(:,1)>100))=0;

        GIAproj = zeros(size(datasets{ii}.datid));
        for jj=1:length(datasets{ii}.siteid)
            sub=find(abs(datasets{ii}.datid- ...
                         datasets{ii}.siteid(jj))<.01);
            GIAproj(sub)=datasets{ii}.siteGIA(jj).*(datasets{ii}.meantime(sub)- ...
                                                    refyear);
        end
    
        datasets{ii}.GIAproj=GIAproj;
        datasets{ii}.Y0=datasets{ii}.Y;
        datasets{ii}.Y=datasets{ii}.Y0-GIAproj;
    else
        datasets{ii}.GIAproj=zeros(size(datasets{ii}.datid));;
        datasets{ii}.Y0=datasets{ii}.Y;
    end
    
end
