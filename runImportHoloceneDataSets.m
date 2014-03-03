% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Mar 2 21:39:29 EST 2014

thinyrs=10;
minlen=50;
firsttime=-1000;

%%%%%%%%%%
datid=[]; time1=[]; time2=[]; mediantime=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[];
siteid=[]; sitenames={}; sitecoords=[];

datPX = importdata(fullfile(IFILES,'RSL_Feb2014.csv'));
datPX.textdata=datPX.textdata(2:end,:);

% catch entries without age errors
sub=find(isnan(datPX.data(:,8))); datPX.data(sub,8)=100;
sub=find(isnan(datPX.data(:,7))); datPX.data(sub,7)=100;



study=datPX.textdata(:,1);
uStudy = unique(study);
for ii=1:length(uStudy)
    sub=find(strcmpi(uStudy{ii},study));
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
    
    
%%%%%%%%%%%%%%%%
    
% Engelhart & Horton database

datHolo=importdata(fullfile(IFILES,'Engelhart_Horton_2012_v4.csv'));
HoloRegions=[1:16];
for curreg=1:length(HoloRegions)
	sub=find((datHolo.data(:,1)==HoloRegions(curreg)).*(datHolo.data(:,2)==0));

	count=[1:length(sub)]';
	wdatid = ones(length(sub),1)*(1e6+HoloRegions(curreg)*1e3);
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

        sitecoords=[sitecoords; mean(wlat) mean(wlong)];

end
siteid=[siteid ; [1e6+HoloRegions*1e3]'];
for i=1:length(HoloRegions)
    sitenames={sitenames{:}, ['EH12_' num2str(HoloRegions(i))]};
    sub=find(datid==(1e6+HoloRegions(i)*1e3));
    sitecoords = [sitecoords; mean(lat(sub)) mean(long(sub))];
end

%%%%

PX.datid=datid;
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
PX.siteid=siteid;
PX.sitenames=sitenames;
PX.meantime=(PX.time1+PX.time2)/2;
PX.sitecoords = sitecoords;

PX0=PX;
sub=find(PX.time1>=firsttime);
shortenfields={'datid','time1','time2','meantime','limiting','Y','dY','compactcorr','istg','lat','long'};
for jj=1:length(shortenfields)
    PX.(shortenfields{jj}) =  PX.(shortenfields{jj})(sub);
end
PX.Ycv=sparse(diag(PX.dY.^2));
PX.Ycv0=sparse(diag(PX.dY.^2));


%%%%%%%%%%%%%%%

% First load tide gauge data

optimizemode=1.0;
[TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords,[],[],[],[],[],optimizemode);

sub=find(TG.datid~=0); subS = find(TG.siteid~=0);
TGNOCW=TG;
shortenfields={'datid','time1','time2','meantime','limiting','Y','dY','compactcorr','istg','lat','long'};
for jj=1:length(shortenfields)
    TGNOCW.(shortenfields{jj}) =  TGNOCW.(shortenfields{jj})(sub);
end
shortenfields={'siteid','sitenames','sitelen'};
for jj=1:length(shortenfields)
    TGNOCW.(shortenfields{jj}) =  TGNOCW.(shortenfields{jj})(subS);
end
TGNOCW.sitecoords=TGNOCW.sitecoords(subS,:);
TGNOCW.Ycv=sparse(TGNOCW.Ycv(sub,sub));
TGNOCW.Ycv0=sparse(TGNOCW.Ycv);

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
[GISfplong,si]=sort(GISfplong); GISfp=GISfp(:,si);
GISfplong=[GISfplong(end)-360 ; GISfplong];
GISfp=[GISfp(:,end) GISfp];
GISfp=GISfp*1000;

%%%%%%

clear datasets;
datasets{1}=MergeDataStructures(TG,PX);
datasets{2}=MergeDataStructures(TGNOCW,PX);
datasets{3}=PX;
datasets{4}=TG;
datasets{5}=TGNOCW;

datasets{1}.label='TG+GSL+PX';
datasets{2}.label='TG+PX';
datasets{3}.label='PX';
datasets{4}.label='TG+GSL';
datasets{5}.label='TG';


for ii=1:length(datasets)
    t1=datasets{ii}.time1; t2=datasets{ii}.time2;
    datasets{ii}.long = mod(datasets{ii}.long,360); sub=find(datasets{ii}.long>180); datasets{ii}.long(sub)=datasets{ii}.long(sub)-360;
    datasets{ii}.meantime=mean([t1 t2],2);
    datasets{ii}.dt = abs(t1-t2)/4;
    datasets{ii}.dY0=datasets{ii}.dY;
    datasets{ii}.Ycv0 = datasets{ii}.Ycv;
    datasets{ii}.compactcorr=sparse(datasets{ii}.compactcorr);


    datasets{ii}.obsGISfp = interp2(GISfplong,GISfplat,GISfp,datasets{ii}.long,datasets{ii}.lat,'linear');
    datasets{ii}.obsGISfp(find(datasets{ii}.lat>100))=1;
    
    % subtract GIA model
    
    ider = datasets{ii}.lat*1e5+datasets{ii}.long;
    [ideru,iderui]=unique(ider);
    
%    [regionsu,regionsusi]=unique(datasets{ii}.datid);
%    sitecoords=[datasets{ii}.lat(regionsusi) datasets{ii}.long(regionsusi)];
%    GIAproju=zeros(size(regionsu));
%    GIAproj=zeros(size(datasets{ii}.Y));
    GIAproju=zeros(size(ideru));
    GIAproj=zeros(size(datasets{ii}.Y));

    for jj=1:length(ideru)
        if datasets{ii}.datid(iderui(jj))>0
            GIAproju(jj)=interp2(ICE5Glat,ICE5Glon,ICE5Ggia,datasets{ii}.lat(iderui(jj)),datasets{ii}.long(iderui(jj)));
            sub=find(ider==ideru(jj));
            GIAproj(sub)=GIAproju(jj).*(datasets{ii}.meantime(sub)-refyear);
        end
    end
    for jj=1:length(datasets{ii}.siteid)
        datasets{ii}.siteGIA(jj)=0;
        sub=find(datasets{ii}.datid == datasets{ii}.siteid(jj));
        if length(sub)>0
            datasets{ii}.siteGIA(jj) = GIAproj(sub(1));
        end
    end
    datasets{ii}.GIAproj=GIAproj;
    datasets{ii}.Y0=datasets{ii}.Y;
    %datasets{ii}.Y=datasets{ii}.Y0-GIAproj;
end



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



