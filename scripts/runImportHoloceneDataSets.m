% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Mar 09 12:48:58 EDT 2015

defval('firsttime',-1300);

thinyrs=10;
minlen=50;
minperstudy=1;
refyear=2010;

%%%%%%%%%%
datid=[]; time1=[]; time2=[]; mediantime=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[];
siteid=[]; sitenames={}; sitecoords=[];

datPX = importdata(fullfile(IFILES,'RSL_Feb2015b.csv'));
datPX.textdata=datPX.textdata(2:end,:);

% catch entries without age errors
sub=find(isnan(datPX.data(:,8))); datPX.data(sub,8)=100;
sub=find(isnan(datPX.data(:,7))); datPX.data(sub,7)=100;

study=datPX.textdata(:,1);
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
            wdY=abs(wY2-wY1)/4;
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

%%%%%%%%%%%%%%%%


datLH = importdata(fullfile(IFILES,'RSL_LateHolocene_Feb2015.csv'));
datLH.textdata=datLH.textdata(2:end,:);

% catch entries without age errors
sub=find(isnan(datLH.data(:,8))); datLH.data(sub,8)=100;
sub=find(isnan(datLH.data(:,7))); datLH.data(sub,7)=100;

study=datLH.textdata(:,1);
uStudy = unique(study);
for ii=1:length(uStudy)
    sub=find(strcmpi(uStudy{ii},study));
    if length(sub)>=minperstudy
        site = datLH.textdata(sub,2);
        uSite=unique(site);
        for jj=1:length(uSite)
            curid = 1e6 + 1e4*ii + jj;
            curstudysite=['LR-' uStudy{ii} '-' uSite{jj}];
            sub2=sub(find(strcmpi(uSite{jj},site)));
            wdatid=ones(length(sub2),1)*curid;
            wmediantime=datLH.data(sub2,6);
            wtime1=wmediantime-datLH.data(sub2,8)+(1:length(sub2))'/1e5;
            wtime2=wmediantime+datLH.data(sub2,7)+(1:length(sub2))'/1e5;
            
            wlimiting=zeros(length(sub2),1);
            wYmedian=datLH.data(sub2,3);
            wY1=wYmedian-datLH.data(sub2,5);
            wY2=wYmedian+datLH.data(sub2,4);
            wY=(wY1+wY2)/2;
            wdY=abs(wY2-wY1)/4;
            wcompactcorr=zeros(length(sub2),1);;
            wlat=datLH.data(sub2,1);
            wlong=datLH.data(sub2,2);
            
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

%%%%

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

% drop too old
sub=find(PX.time1>=firsttime);
sub=intersect(sub,find(abs(PX.lat)<=90));
subS=find(abs(PX.sitecoords(:,1))<=90);

PX=SubsetDataStructure(PX,sub,subS);


sub=find(PX.datid<1e6);
subS=find(PX.siteid<1e6);
PXnoLH=SubsetDataStructure(PX,sub,subS);

sub=find(PX.datid>=1e6);
subS=find(PX.siteid>=1e6);
LH=SubsetDataStructure(PX,sub,subS);

%%%%%%%
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
        dat=importdata(fullfile(IFILES,'amsterdam.sea.level.txt'));
        wsite='AMSTERDAM_OLD'; wcurid=5001;
        wlat = [52.3667]; wlong=[4.9000];
        wdatid = ones(length(dat.data),1)*wcurid;
        wtime = dat.data(:,1);
        wY = dat.data(:,2);
    elseif ppp==2
        dat=importdata(fullfile(IFILES,'Kronstadt_ReportsFGI_Bogdanov_appendix.csv'));
        wsite='KRONSTADT_OLD'; wcurid=5002;
        wlat = [59.98]; wlong=[29.77];
        wdatid = ones(length(dat.data),1)*wcurid;
        wtime = dat.data(:,1);
        wY = dat.data(:,3);
    elseif ppp==3
         dat=importdata(fullfile(IFILES,'ekman_2003_stockholm.csv'));
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

TGoldS = GPSmoothTideGauges(TGold,11,1.0,10,[],1700);
%%%%%%%%%%%%%%%

% First load tide gauge data
optimizemode=1.0;

psmsldir=fullfile(IFILES,'rlr_annual');

[TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords,[],[],[],[],[],optimizemode,psmsldir,'none');

% account for additional uncertainties in GSL curve;

gslfile=fullfile(IFILES,'Hay2014_KFandGP_GMSL.mat');

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
TGNOGSL=TG;
TG=MergeDataStructures(TG,HayGSL);


%%%%%%%%%%%%%%%%



%%%%%% load GIA model

giamodel.gia=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Dsea_250');
giamodel.lat=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lat');
giamodel.long=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lon');

ICE5Ggia=giamodel.gia;
ICE5Glat=giamodel.lat;
ICE5Glon=giamodel.long;
sub=find(ICE5Glon>180); ICE5Glon(sub)=ICE5Glon(sub)-360;
[ICE5Glon,si]=sort(ICE5Glon); ICE5Ggia=ICE5Ggia(si,:);

% fingerprint

[GISfp,GISfplong,GISfplat]=readFingerprintInd('gis',IFILES);
sub=find(GISfplong>180); GISfplong(sub)=GISfplong(sub)-360;
[GISfplong,si]=sort(GISfplong); GISGISfp=GISfp(:,si);
GISfplong=[GISfplong(end)-360 ; GISfplong];
GISfp=[GISfp(:,end) GISfp];
GISfp=GISfp*1000;
%%%%%%%

dat=importdata(fullfile(IFILES,'Grinsted2009_JonesA1B.txt'));
Grinsted2009_Jones.year=dat.data(:,1);
Grinsted2009_Jones.quantiles=[5 16 50 84 95];
Grinsted2009_Jones.y=dat.data(:,2:end)*1000;


dat=importdata(fullfile(IFILES,'Grinsted2009_MobergA1B.txt'));
Grinsted2009_Moberg.year=dat.data(:,1);
Grinsted2009_Moberg.quantiles=[5 16 50 84 95];
Grinsted2009_Moberg.y=dat.data(:,2:end)*1000;

sub=find(mod(Grinsted2009_Moberg.year,20)==0);
sub=intersect(sub,find(Grinsted2009_Moberg.year<=2010));
Grinsted.Y=Grinsted2009_Moberg.y(sub,3);
Grinsted.dY=(Grinsted2009_Moberg.y(sub,5)-Grinsted2009_Moberg.y(sub,1))/norminv(.95)/2;
Grinsted.Ycv=(diag(Grinsted.dY)).^2;
Grinsted.datid=0*Grinsted.Y;
Grinsted.time1=Grinsted2009_Moberg.year(sub);
Grinsted.time2=Grinsted.time1;
Grinsted.meantime=Grinsted.time1;
Grinsted.lat=ones(size(Grinsted.Y))*1e6;
Grinsted.long=Grinsted.lat;
Grinsted.compactcorr=0*Grinsted.Y;
Grinsted.limiting=0*Grinsted.Y;
Grinsted.istg=ones(size(Grinsted.Y));
Grinsted.siteid=0;
Grinsted.sitenames={'Grinsted_Moberg'};
Grinsted.sitecoords=[1e6 1e6];
Grinsted.sitelen=length(Grinsted.Y);

%%%%%%

clear datasets;
datasets{1}=MergeDataStructures(TG,PX);
datasets{2}=MergeDataStructures(MergeDataStructures(TG,PX),GSLflattener);
datasets{3}=MergeDataStructures(TG,PXnoLH);
datasets{4}=PX;

datasets{1}.label='TG+GSL+PX';
datasets{2}.label='TG+GSL+PX+flat';
datasets{3}.label='TG+PXnoLH';
datasets{4}.label='PX';


for ii=1:length(datasets)
    t1=datasets{ii}.time1; t2=datasets{ii}.time2;
    datasets{ii}.long = mod(datasets{ii}.long,360); sub=find(datasets{ii}.long>180); datasets{ii}.long(sub)=datasets{ii}.long(sub)-360;
    datasets{ii}.meantime=mean([t1 t2],2);
    datasets{ii}.dt = abs(t1-t2)/4;
    datasets{ii}.compactcorr=sparse(datasets{ii}.compactcorr);


    datasets{ii}.obsGISfp = interp2(GISfplong,GISfplat,GISfp,datasets{ii}.long,datasets{ii}.lat,'linear');
    datasets{ii}.obsGISfp(find(datasets{ii}.lat>100))=1;
    datasets{ii}.siteGISfp = interp2(GISfplong,GISfplat,GISfp,datasets{ii}.sitecoords(:,2),datasets{ii}.sitecoords(:,1),'linear');
    datasets{ii}.siteGISfp(find(datasets{ii}.sitecoords(:,1)>100))=1;
    
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
end



% indices

datNAO = importdata(fullfile(IFILES,'nao-trouet2009.txt'));

%% comparison to temperature records
impt=importdata(fullfile(IFILES,'Marcott2013_global.txt'));
Marcottgl.yr = 1950-impt.data(:,1);
Marcottgl.T = impt.data(:,2);
Marcottgl.dT = impt.data(:,3);

impt=importdata(fullfile(IFILES,'Marcott2013_regional.txt'));
Marcottr.yr=1950-impt.data(:,1);
Marcott.Next_T=impt.data(:,2);
Marcott.Next_dT=impt.data(:,3);
Marcott.eq_T=impt.data(:,4);
Marcott.eq_dT=impt.data(:,5);
Marcott.Sext_T=impt.data(:,6);
Marcott.Sext_dT=impt.data(:,7);

impt=importdata(fullfile(IFILES,'HadCRUT4-gl.dat'));
HadCRUT.yr=impt(1:2:end,1);
HadCRUT.T=impt(1:2:end,end);

Mann.T=importdata(fullfile(IFILES,'glglfulihad_smxx.txt'));
Mann.yr = [[1:length(Mann.T)]-1]';

impt=importdata(fullfile(IFILES,'Lund2006_transport.csv'));
Lundtransport.yr=impt.data(:,1);
Lundtransport.Sv=impt.data(:,2);
Lundtransport.dSv=impt.data(:,3);


impt=importdata(fullfile(IFILES,'Cronin2010_ChesapeakeT.txt'));
Chesapeake.yr = impt(:,1);
Chesapeake.T = impt(:,2);

impt=importdata(fullfile(IFILES,'Richey2007_GulfOfMexico.csv'));
GOM.yr = 1950-impt.data(:,1);
GOM.T = impt.data(:,4);

impt=importdata(fullfile(IFILES,'RothJoos2013TSI.csv'));
RothJoos.yr = impt.data(:,1);
RothJoos.TSI = impt.data(:,2);
RothJoos.TSIerr = impt.data(:,3);

impt=importdata(fullfile(IFILES,'haug2001_cariaco_ti.txt'));
haug.yr = 1950-impt.data(:,1);
haug.Ti = impt.data(:,2);


dat=importdata(fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv'));
CW2011.year=dat.data(:,1);
CW2011.y=dat.data(:,2);
CW2011.dy=dat.data(:,3);
