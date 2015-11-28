% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Nov 27 13:07:50 EST 2015

wdataset=datasets{2};

%firstyears=[0 0 700 1400 1800 1900];
%lastyears=[1700 700 1400 1800 1900 2000];

firstyears=[0];
lastyears=[2000];

for www=1:length(firstyears)

    clf;
    ax = worldmap('World');
    setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off' );
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;

    ud=unique(wdataset.datid(find((wdataset.meantime>firstyears(www)).*(wdataset.meantime<lastyears(www)))));
    ud=ud(ud>0);
    sub=find(ismember(wdataset.siteid,ud));
    subtg=intersect(sub,find(ismember(wdataset.siteid,wdataset.datid(find(wdataset.istg)))));
    subpx=setdiff(sub,subtg);
    clear hs1 hs2;
    %hs2=scatterm(wdataset.sitecoords(subtg,1),wdataset.sitecoords(subtg,2),15,'g','filled','MarkerFaceColor','g','Marker','s','MarkerEdgeColor','g'); hold on;
 hs1=scatterm(wdataset.sitecoords(subpx,1),wdataset.sitecoords(subpx,2),15,'r','filled','MarkerFaceColor','r','Marker','d','MarkerEdgeColor','r'); hold on;
        
    axis tight;
    %title([num2str(firstyears(www)) '-' num2str(lastyears(www))]);
    ht=title('a) Proxy data');
    set(ht,'fontsize',12);
    pdfwrite(['sitemap_proxy_' num2str(firstyears(www)) '_' num2str(lastyears(www))]);


    clf;
    ax = worldmap('World');
    setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off' );
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
    hold on;
   hs2=scatterm(wdataset.sitecoords(subtg,1),wdataset.sitecoords(subtg,2),15,'g','filled','MarkerFaceColor','g','Marker','s','MarkerEdgeColor','g'); hold on;
   %hs1=scatterm(wdataset.sitecoords(subpx,1),wdataset.sitecoords(subpx,2),15,'b','filled','MarkerFaceColor','b','Marker','d','MarkerEdgeColor','b'); hold on;
        
    axis tight;
    %    title([num2str(firstyears(www)) '-' num2str(lastyears(www))]);

       %title([num2str(firstyears(www)) '-' num2str(lastyears(www))]);
    ht=title('b) Tide-gauge data');
    set(ht,'fontsize',12);
 pdfwrite(['sitemap_tidegauge_' num2str(firstyears(www)) '_' num2str(lastyears(www))]);

end
