% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Nov 9 2014

defval('firsttime',-1000);

thinyrs=10;
minlen=50;
minperstudy=2;

%%%%%%%%%%
datid=[]; time1=[]; time2=[]; mediantime=[]; limiting=[]; Y=[]; dY = []; compactcorr = [];
istg = []; lat=[]; long=[];
siteid=[]; sitenames={}; sitecoords=[];

datPX = importdata(fullfile(IFILES,'RSL_Nov2014b.csv'));
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
    wcompactcorr = wY;
    
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
    sitenames={sitenames{:}, ['EH12_' num2str(HoloRegions(curreg))]};
    siteid=[siteid ; [1e6+curreg*1e3]'];
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

% drop too old, and too near field
sub=find(PX.time1>=firsttime);
sub=intersect(sub,find(abs(PX.lat)<=80));
subS=find(abs(PX.sitecoords(:,1))<=80);

PX=SubsetDataStructure(PX,sub,subS);

sub = find(PX.datid<1e6);
subS = find(PX.siteid<1e6);
PXnoEH = SubsetDataStructure(PX,sub,subS);

% drop Greenland and Iceland

subexcl = union(strmatch('Greenland',PX.sitenames),strmatch('Iceland',PX.sitenames));
sub = find(~ismember(PX.datid,PX.siteid(subexcl)));
subS = find(~ismember(PX.siteid,PX.siteid(subexcl)));
PXnonf = SubsetDataStructure(PX,sub,subS);

% $$$ incls = {{'North Carolina','New Jersey','Florida'},{'Florida','New Jersey','North Carolina','Nova Scotia'},{'Florida','New Jersey','North Carolina','Nova Scotia','New Zealand','Cook Islands'},{'Florida','New Jersey','North Carolina','Nova Scotia','New Zealand','Cook Islands','Massachusetts','Maine','Louisiana'}};
% $$$ 
% $$$ for ii=1:length(incls);
% $$$     incl=incls{ii};
% $$$     subincl=[]; subinclS=[];
% $$$     for jj=1:length(incl)
% $$$         doincl=find(strncmp(incl{jj},PX.sitenames, ...
% $$$                             length(incl{jj})));
% $$$         subinclS=union(subinclS,doincl);
% $$$         for kk=1:length(doincl)
% $$$             subincl=union(subincl,find(PX.datid==PX.siteid(doincl(kk))));
% $$$         end
% $$$     end
% $$$     sub=subincl; subS=subinclS;
% $$$     PXsub{ii}=SubsetDataStructure(PX,sub,subS);
% $$$ end

% now create a slimmed version for training
% $$$ excl={'France','Tasmania','Spain','New Zealand','Italy','Israel',['Isle ' ...
% $$$                     'of Wight'],'Iceland','Greenland'};
% $$$ 
% $$$ subexcl=[]; subexclS=[];
% $$$ 
% $$$ for ii=1:length(excl)
% $$$     doexcl=find(strncmp(excl{ii},PX.sitenames,length(excl{ii})));
% $$$     subexclS=union(subexclS,doexcl);
% $$$     for jj=1:length(doexcl)
% $$$         subexcl=union(subexcl,find(PX.datid==PX.siteid(doexcl(jj))));
% $$$     end
% $$$ end
% $$$ sub=setdiff(1:length(PX.datid),subexcl);
% $$$ subS=setdiff(1:length(PX.siteid),subexclS);
% $$$ 
% $$$ PXslim=SubsetDataStructure(PX,sub,subS);


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
gslfile=fullfile(IFILES,'/CSIRO_Recons_gmsl_yr_2011.csv');

[TG,TG0,thetL,TGmodellocal] = GPSmoothNearbyTideGauges(PX.sitecoords,[],[],[],[],[],optimizemode,psmsldir,gslfile);

% account for additional uncertainties in GSL curve;

sub=find(TG.datid==0);
GSLnewcv = TG.Ycv(sub,sub);

% reduced degrees of freedom (1/4 year)
GSLnewcv = GSLnewcv * 4;

% 2s trend error of 0.1 mm/y over a century due to GIA
GSLtrenderror = 0.1/2;
GSLnewcv = GSLnewcv + bsxfun(@times,(TG.meantime(sub)-2000),(TG.meantime(sub)'-2000)) * (GSLtrenderror)^2;
TG.dY(sub)=sqrt(diag(GSLnewcv));
TG.Ycv(sub,sub)=GSLnewcv;


% drop near field
sub=union(find(abs(TG.lat)<=53),find(TG.lat>100));
subS=union(find(abs(TG.sitecoords(:,1))<=53),find(TG.sitecoords(:,1)>100));

TG=SubsetDataStructure(TG,sub,subS);

% add in old tide gauges
TG=MergeDataStructures(TG,TGoldS);

%

sub=find(TG.datid~=0); subS = find(TG.siteid~=0);
TGNOCW=SubsetDataStructure(TG,sub,subS);


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

%%%%%%

clear datasets;
datasets{1}=MergeDataStructures(TGNOCW,PX);
datasets{2}=PX;
datasets{3}=MergeDataStructures(TGNOCW,PXnoEH);
datasets{4}=MergeDataStructures(TG,PXnoEH);
datasets{5}=MergeDataStructures(TG,PX);
datasets{6}=MergeDataStructures(TG,PXnonf);
datasets{7}=MergeDataStructures(TGNOCW,PXnonf);

datasets{1}.label='TG+PX';
datasets{2}.label='PX';
datasets{3}.label='TG+PXnoEH';
datasets{4}.label='TG+GSL+PXnoEH';
datasets{5}.label='TG+GSL+PX';
datasets{6}.label='TG+GSL+PXnonf';
datasets{7}.label='TG+PXnonf';

%for jj=1:length(PXsub)
%    datasets{end+1}=MergeDataStructures(TGNOCW,PXsub{jj});
%    datasets{end}.label=['TG+PXsub' num2str(jj)];
%end


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
