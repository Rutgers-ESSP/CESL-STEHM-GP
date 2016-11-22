% Read in data files
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Nov 21 21:59:10 EST 2016

minperstudy=1; % minimum size of studies to include

loadNearbyTG=1; % load nearby tide gauge sites?
psmsldir=fullfile(IFILES,'rlr_annual'); % location of PSMSL data files
thinyrs=10; % how much to thin/smooth tide gauge data
optimizemode=1.0; % local optimization only for smoothing GPs for tide gauges

usePriorGIA=1; % use a GIA model as the mean of the prior?
refyear=2010; % reference year for GIA hindcast

%%%%%%%%%%
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

% if you want to subset in any way, do it here
% currently just subsetting those with latitude <= 90 (i.e., everything)
sub=find(abs(PX.lat)<=90);
subS=find(abs(PX.sitecoords(:,1))<=90);

PX=SubsetDataStructure(PX,sub,subS);

%%%%%%%

% Load nearby tide gauge data
TGnoisemask = [];
if loadNearbyTG
    [TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords,[],[],[],[],[],optimizemode,psmsldir,'none',[],TGnoisemask);
end
    
%%%%%% load GIA model

if usePriorGIA
    giamodel.gia=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Dsea_250');
    giamodel.lat=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lat');
    giamodel.long=ncread(fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'),'Lon');

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

%%%%%%%

% create data structures

clear datasets;
datasets{1}=PX;
datasets{1}.label='PX';

if loadNearbyTG
    datasets{2}=MergeDataStructures(TG,PX);
    datasets{2}.label='TG+PX';
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