% Generate maps of rate fields.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2017-07-25 18:18:44 -0400

% null data set for prior cvalues
nulldataset=SubsetDataStructure(wdataset,1,1);
nulldataset.meantime=2000; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;

trainsub = find((wdataset.limiting==0));
scattersize=300;
latlim=[25 50];
longlim=[275 305];
gridsz=[.25 .25];

firstyears2=[ 0     0   400 800   1200 1600 1800 1900];
lastyears2= [ 1700  400 800 1200  1600 1800 1900 2000];
doNoiseMasks=ones(size(firstyears2));
cutoffThreshold=ones(size(firstyears2))*sqrt(.8); 

% initialize grid

Flat=latlim(1):gridsz(1):latlim(2);
Flong=longlim(1):gridsz(2):longlim(2);
sub=find(testsites(:,2)<=360);
ulat=unique(round(testsites(sub,2)));
ulong=unique(round(testsites(sub,3)));
ulat=unique(bsxfun(@plus,ulat,[-4:gridsz:4]));
ulong=unique(mod(bsxfun(@plus,ulong,[-4:gridsz:4]),360));
sub=find(abs(ulat)<90); ulat=ulat(sub); ulat=ulat(:)'; ulong=ulong(:)';
Flat=union(Flat,ulat);
Flong=union(mod(Flong,360),mod(ulong,360));

Flat1=min(Flat):max(Flat);
Flong1=min(Flong):max(Flong);

[FLONG,FLAT]=meshgrid(Flong,Flat);
[FLONG1,FLAT1]=meshgrid(Flong1,Flat1);

land = shaperead('landareas', 'UseGeoCoords', true);

covered=[];
for sss=1:length(land)
    covered=union(covered,find(inpolygon(FLAT1(:),FLONG1(:),land(sss).Lat,land(sss).Lon)));
    %disp(sprintf('%0.0f: %0.0f',[sss length(covered)]));     
end
uncovered=setdiff(1:length(FLAT1(:)),covered);

GSLsitesub=find(testsites==0);
GSLdatsub=find(testreg==0);

figure;

% %%%%

for qqq=1:length(firstyears2)
    firstyears=firstyears2(qqq);
    lastyears=lastyears2(qqq);
    doNoiseMask=doNoiseMasks(qqq);

    [fslopeF,sdslopeF,~,~,~,~,~,~,passderivs,invcv] = RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(doNoiseMask,:),Flat,Flong,firstyears,lastyears,trainsub,ICE5G);    
    [fslopeGSL,sdslopeGSL]=SLRateCompare(f2s{iii}(GSLdatsub,1),V2s{iii}(GSLdatsub,GSLdatsub,1),testsites(GSLsitesub),testreg(GSLdatsub),testX(GSLdatsub,3),firstyears,lastyears);
    [priorslope,sdpriorslope] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,0,firstyears,lastyears);    

    mapped = griddata(FLONG(:),FLAT(:),fslopeF,Flong1,Flat1(:),'linear');
    sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF,Flong1,Flat1(:),'linear');

    u=sdmapped/sdpriorslope;
    threshold=cutoffThreshold(qqq);
    subbad=find((sdmapped/sdpriorslope)>threshold);
    mapped(subbad)=NaN;

    %%

    clf;
    ax = worldmap(latlim,longlim);
    setm(ax, 'meridianlabel','off','parallellabel','off','flinewidth',3);
    %geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;

    subgood=find(~isnan(mapped));
    subgood=intersect(subgood,uncovered);
    hs1=scatterm(FLAT1(subgood),FLONG1(subgood),scattersize,mapped(subgood),'filled','marker','s');

    hold on;

    sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
    subbad2=intersect(subbad,sublong);

    colormap(parula);
    axis tight;
    hcb=colorbar;
    set(hcb,'fontsize',12);

    box on;
    tl=get(hcb,'ticks');
    tlab=get(hcb,'ticklabels');
    for sss=1:length(tl)
        if abs((tl(sss)*10)-round(tl(sss)*10))>.01
            tlab{sss}='';
        else
            tlab{sss}=sprintf(' %0.1f',tl(sss));
        end
        
    end
    set(hcb,'ticklabels',tlab);

    umap=geoshow(ax, land, 'FaceColor',  [0.85 0.85 0.85]);

    ht=title([num2str(firstyears2(qqq)) '-' num2str(lastyears2(qqq)) ' CE']);
    %set(ht,'fontsize',16);
    pdfwrite(['fieldmap_' labl '_' num2str(firstyears2(qqq)) '_' num2str(lastyears2(qqq)) '_' noiseMasklabels{doNoiseMask}]);

    %%% now plot standard deviation

    clf;
    ax = worldmap(latlim,longlim);
    setm(ax, 'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;

    hs1=scatterm(FLAT1(:),FLONG1(:),scattersize,sdmapped(:),'filled','marker','s');

    hold on;
    colormap(parula);
    axis tight;
    hcb=colorbar;

    box on;
    %   caxis([0 1.5]);
    geoshow(ax, land, 'FaceColor', 'none');
    title([num2str(firstyears2(qqq)) '-' num2str(lastyears2(qqq)) ' (mm/y)']);
    pdfwrite(['fieldmap_' labl '_' num2str(firstyears2(qqq)) '_' num2str(lastyears2(qqq)) '_sd' '_' noiseMasklabels{doNoiseMask}]);

end