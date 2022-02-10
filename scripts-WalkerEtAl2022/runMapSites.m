% Generate maps of proxy and tide gauge data.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:07:40 -0400

wdataset=datasets{2};

firstyears=[0];
lastyears=[2000];
www=1;

% proxy data

clf;
ax = worldmap(latlim,longlim);
setm(ax, 'meridianlabel','off','parallellabel','off' );
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on;

ud=unique(wdataset.datid(find((wdataset.meantime>firstyears(www)).*(wdataset.meantime<lastyears(www)))));
ud=ud(ud>0);
sub=find(ismember(wdataset.siteid,ud));
subtg=intersect(sub,find(ismember(wdataset.siteid,wdataset.datid(find(wdataset.istg)))));
subpx=setdiff(sub,subtg);
clear hs1 hs2;
hs1=scatterm(wdataset.sitecoords(subpx,1),wdataset.sitecoords(subpx,2),15,'r','filled','MarkerFaceColor','r','Marker','d','MarkerEdgeColor','r'); hold on;

axis tight;
ht=title('a) Proxy data');
set(ht,'fontsize',12);
pdfwrite(['sitemap_proxy_' num2str(firstyears(www)) '_' num2str(lastyears(www))]);


% tide gauge data

clf;
ax = worldmap(latlim,longlim);
setm(ax,'meridianlabel','off','parallellabel','off' );
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on;
hs2=scatterm(wdataset.sitecoords(subtg,1),wdataset.sitecoords(subtg,2),15,'g','filled','MarkerFaceColor','g','Marker','s','MarkerEdgeColor','g'); hold on;

axis tight;
ht=title('b) Tide-gauge data');
set(ht,'fontsize',12);
pdfwrite(['sitemap_tidegauge_' num2str(firstyears(www)) '_' num2str(lastyears(www))]);
