% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 14:50:18 EST 2014

   nulldataset=SubsetDataStructure(wdataset,1,1);
    nulldataset.meantime=2000; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;
    
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
    
    [fslopeF,sdslopeF,~,~,~,~,~,~,passderivs,invcv] = RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears,lastyears,trainsub,ICE5G,passderivs,invcv);    
    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);
    [priorslope,sdpriorslope] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,0,firstyears,lastyears);    
  
    clf;
    ax = worldmap('World');
    setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off' );
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
    
    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);
    
    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    mapped = griddata(FLONG(:),FLAT(:),fslopeF,Flong1,Flat1(:),'linear');
    sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF,Flong1,Flat1(:),'linear');
    
    u=sdmapped/sdpriorslope;
    threshold=max(.35,quantile(u(:),.05));
    subbad=find((sdmapped/sdpriorslope)>threshold);
    mapped(subbad)=NaN;
    
    hs1=scatterm(FLAT1(:),FLONG1(:),10,mapped(:),'filled','marker','s');
    %hs2=plotm(FLAT1(subbad),FLONG1(subbad),'o','color',[1 1 1],'markersize',4,'markerfacecolor','w','markeredgecolor','w')
    
    hold on;
    geoshow(ax, land, 'FaceColor', 'none');

    ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(qqq)).*(wdataset.time1<=lastyears(qqq)))));
    sub1=find(ismember(wdataset.siteid,ud));
    %hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
    
    axis tight;
    hcb=colorbar;

    box on;
    %   caxis([0 1.5]);
    caxis([-1 3]);
    title({[num2str(firstyears(qqq)) '-' num2str(lastyears(qqq)) ' (mm/y)'],['GSL: ' sprintf('%0.2f \\pm %0.2f mm/y',[fslopeavg(1) 2*sdslopeavg(1)])  ' --- threshold \sigma = ' sprintf('%0.2f',threshold*sdpriorslope) ' mm/y' ]});
    pdfwrite(['fieldmap_' labl '_' num2str(firstyears(qqq)) '_' num2str(lastyears(qqq))]);
    
    %%%%
    
    firstyears2=[-500 1000 1500 1800 1900];
    lastyears2=[1000 1500 1800 1900 2000];
    
    for qqq=1:length(firstyears2)
        firstyears=[0 firstyears2(qqq)  ];
        lastyears=[1800 lastyears2(qqq)];    
        
        [fslopeF,sdslopeF,~,~,fslopediffF,sdslopediffF,diffplusF,difflessF,passderivs,invcv] = RegressRateField(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears,lastyears,trainsub,ICE5G,passderivs,invcv);    
     [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);
         [~,~,~,~,priorslope,sdpriorslope] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,0,firstyears,lastyears);    
 
       
        clf;
        ax = worldmap('World');
        setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off' );
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
        hold on;
        
        Flat1=min(Flat):max(Flat);
        Flong1=min(Flong):max(Flong);
        
        [FLONG,FLAT]=meshgrid(Flong,Flat);
        [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
        mapped = griddata(FLONG(:),FLAT(:),fslopediffF,Flong1,Flat1(:),'linear');
        sdmapped = griddata(FLONG(:),FLAT(:),sdslopediffF,Flong1,Flat1(:),'linear');
        
        u=sdmapped/sdpriorslope;
        threshold=max(.35,quantile(u(:),.05));
        subbad=find((sdmapped/sdpriorslope)>threshold);
        mapped(subbad)=NaN;
        
        hs1=scatterm(FLAT1(:),FLONG1(:),10,mapped(:),'filled','marker','s');
        %hs2=plotm(FLAT1(subbad),FLONG1(subbad),'o','color',[1 1 1],'markersize',4,'markerfacecolor','w','markeredgecolor','w')
        hold on;
        geoshow(ax, land, 'FaceColor', 'none');

        ud=unique(wdataset.datid(find((wdataset.time2>=firstyears2(qqq)).*(wdataset.time1<=lastyears2(qqq)))));
        sub1=find(ismember(wdataset.siteid,ud));
        %hs1=scatterm(wdataset.sitecoords(sub1,1),wdataset.sitecoords(sub1,2),15,'k','filled','MarkerFaceColor','w','Marker','d','MarkerEdgeColor','k'); hold on;
        
        axis tight;
        hcb=colorbar;

        box on;
        %   caxis([0 1.5]);
        q=caxis;
        title({[num2str(firstyears2(qqq)) '-' num2str(lastyears2(qqq)) ' relative to 0-1800 (mm/y)'],['GSL: ' sprintf('%0.2f \\pm %0.2f mm/y',[fslopeavgdiff(1) 2*sdslopeavgdiff(1)]) ' --- threshold \sigma = ' sprintf('%0.2f',threshold*sdpriorslope) ' mm/y' ]});
        caxis([min(0,q(1)) q(2)]);
        
        pdfwrite(['fieldmap_rel01800_' labl '_' num2str(firstyears2(qqq)) '_' num2str(lastyears2(qqq))]);
    end