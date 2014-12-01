% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 30 22:12:33 EST 2014

nulldataset=SubsetDataStructure(wdataset,1,1);
nulldataset.meantime=2000; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;

trainsub = find((wdataset.limiting==0));

%    firstyears=[0   -500 1000 1500 1800 1900];
%lastyears=[1800 1000 1500 1800 1900 2000];
qqq=1;
firstyears=[0  ];
lastyears=[1800];
Flat=-85:10:85;
Flong=0:20:360;
sub=find(testsites(:,2)<=360);
ulat=unique(round(testsites(sub,2)));
ulong=unique(round(testsites(sub,3)));
ulat=unique(bsxfun(@plus,ulat,[-4:2:4]));
ulong=unique(mod(bsxfun(@plus,ulong,[-4:2:4]),360));
sub=find(abs(ulat)<90); ulat=ulat(sub); ulat=ulat(:)'; ulong=ulong(:)';
Flat=union(Flat,ulat);
Flong=union(mod(Flong,360),mod(ulong,360));

GSLsitesub=find(testsites==0);
GSLdatsub=find(testreg==0);

[fslopeF,sdslopeF,~,~,~,~,~,~,passderivs,invcv] = RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears,lastyears,trainsub,ICE5G);    
[fslopeGSL,sdslopeGSL]=SLRateCompare(f2s{iii}(GSLdatsub,1),V2s{iii}(GSLdatsub,GSLdatsub,1),testsites(GSLsitesub),testreg(GSLdatsub),testX(GSLdatsub,3),firstyears,lastyears);
[priorslope,sdpriorslope] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,0,firstyears,lastyears);    


clf;
ax = worldmap('World');
setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off','flinewidth',3);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on;

Flat1=min(Flat):max(Flat);
Flong1=min(Flong):max(Flong);

[FLONG,FLAT]=meshgrid(Flong,Flat);
[FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
mapped = griddata(FLONG(:),FLAT(:),fslopeF,Flong1,Flat1(:),'linear')+fslopeGSL;
sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF,Flong1,Flat1(:),'linear');

u=sdmapped/sdpriorslope;
%threshold=max(.35,quantile(u(:),.05));
threshold=sqrt(.8);
subbad=find((sdmapped/sdpriorslope)>threshold);
mapped(subbad)=NaN;

hs1=scatterm(FLAT1(:),FLONG1(:),10,mapped(:),'filled','marker','s');
%hs2=plotm(FLAT1(subbad),FLONG1(subbad),'o','color',[1 1 1],'markersize',4,'markerfacecolor','w','markeredgecolor','w')

hold on;

sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
subbad2=intersect(subbad,sublong);
%hbad=plotm(FLAT1(subbad2),FLONG1(subbad2),'kx');
%set(hbad,'MarkerSize',3,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
%hbad=scatterm(FLAT1(subbad2),FLONG1(subbad2),10,[.5 .5 .5],'x');

geoshow(ax, land, 'FaceColor',  [0.85 0.85 0.85]);

ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(qqq)).*(wdataset.time1<=lastyears(qqq)))));
sub1=find(ismember(wdataset.siteid,ud));
%hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;

colormap(jet);
axis tight;
hcb=colorbar;

box on;
%   caxis([0 1.5]);
caxis([-1 3]);
title({[num2str(firstyears(qqq)) '-' num2str(lastyears(qqq)) ' (mm/y)'],['GSL: ' sprintf('%0.2f \\pm %0.2f mm/y',[fslopeGSL 2*sdslopeGSL])  ' --- threshold \sigma = ' sprintf('%0.2f',threshold) '\sigma_0' ]});
pdfwrite(['fieldmap_' labl '_' num2str(firstyears(qqq)) '_' num2str(lastyears(qqq))]);

%%%%

firstyears2=[ 0   400 800  1200 1600 1800 1900];
lastyears2= [ 400 800 1200 1600 1800 1900 2000]; 
for qqq=1:length(firstyears2)

    firstyears=[0 firstyears2(qqq)];
    lastyears = [1800 lastyears2(qqq)];
    disp(sprintf('%0.0f--%0.0f',[firstyears2(qqq) lastyears2(qqq)]));
    
    [~,~,~,~,fslopediffF,sdslopediffF,diffplusF,difflessF,passderivs,invcv] = RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears,lastyears,trainsub,ICE5G,passderivs,invcv);    
    [~,~,fslopeGSL,sdslopeGSL]=SLRateCompare(f2s{iii}(GSLdatsub,1),V2s{iii}(GSLdatsub,GSLdatsub,1),testsites(GSLsitesub),testreg(GSLdatsub),testX(GSLdatsub,3),firstyears,lastyears);
    [~,~,~,~,priorslope,sdpriorslope] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,0,firstyears,lastyears); 
    %fslopediffF=fslopediffF+fslopeGSL;

    
    clf;
    ax = worldmap('World');
    setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
    
    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);
    
    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    mapped = griddata(FLONG(:),FLAT(:),fslopediffF,Flong1,Flat1(:),'linear');%+fslopeGSL;
    sdmapped = griddata(FLONG(:),FLAT(:),sdslopediffF,Flong1,Flat1(:),'linear');
    
    u=sdmapped/sdpriorslope;
    %    threshold=max(.35,quantile(u(:),.05));
    %threshold=sqrt(.5);
    %subbad=find(u>threshold);
    %mapped(subbad)=NaN;
    
    subbad=find(abs(normcdf(mapped./sdmapped)-.5)<.166);
    
    
    
    cmap=colormap(jet);
    cmap(:,1)=.5*cmap(:,2)+cmap(:,1); cmap(:,3)=.5*cmap(:,2)+cmap(:,3);
    cmap=min(1,cmap); colormap(cmap);
    
    hs1=scatterm(FLAT1(:),FLONG1(:),10,mapped(:),'filled','marker','s');
    %hs2=plotm(FLAT1(subbad),FLONG1(subbad),'o','color',[1 1 1],'markersize',4,'markerfacecolor','w','markeredgecolor','w')
    hold on;
    
    
    % subbad2=find(abs((mapped-fslopeGSL)./sdmapped)<norminv(.667));
    sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));
    subbad2=intersect(subbad,sublong);
    %hbad=plotm(FLAT1(subbad2),FLONG1(subbad2),'kx');
    %set(hbad,'MarkerSize',3,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
    hbad=scatterm(FLAT1(subbad2),FLONG1(subbad2),10,[.5 .5 .5],'x');
    
    geoshow(ax, land, 'FaceColor',  [0.85 0.85 0.85]);

    ud=unique(wdataset.datid(find((wdataset.time2>=firstyears2(qqq)).*(wdataset.time1<=lastyears2(qqq)))));
    sub1=find(ismember(wdataset.siteid,ud));
    %hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
    
    axis tight;
    hcb=colorbar;
    

    box on;
    %   caxis([0 1.5]);
    q1=min(fslopediffF(:)-fslopeGSL);
    q2=max(fslopediffF(:)-fslopeGSL);
    caxis(fslopeGSL+[-1 1]*max(abs(q1),abs(q2)));
    title({[num2str(firstyears2(qqq)) '-' num2str(lastyears2(qqq)) ' relative to 0-1800 (mm/y)'],['GSL: ' sprintf('%0.2f \\pm %0.2f mm/y',[fslopeGSL(1) 2*sdslopeGSL(1)]) ' --- threshold \sigma = ' sprintf('%0.2f',threshold) '\sigma_0' ]});
    
    pdfwrite(['fieldmap_rel01800_' labl '_' num2str(firstyears2(qqq)) '_' num2str(lastyears2(qqq))]);
end
colormap(jet);