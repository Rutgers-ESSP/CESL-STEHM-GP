% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Nov 15 13:40:51 EST 2014

clf;
ax = worldmap('World');
setm(ax, 'Origin',[0 -90 0],'meridianlabel','off','parallellabel','off' );
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on;

ud=unique(wdataset.datid(find((wdataset.time2>=0))));
ud=ud(find(ud>0));
sub=find(ismember(wdataset.siteid,ud));
subtg=intersect(sub,find(ismember(wdataset.siteid,wdataset.datid(find(wdataset.istg)))));
subpx=setdiff(sub,subtg);
clear hs1 hs2;
hs1=scatterm(wdataset.sitecoords(subpx,1),wdataset.sitecoords(subpx,2),15,'k','filled','MarkerFaceColor','b','Marker','d','MarkerEdgeColor','k'); hold on;
hs2=scatterm(wdataset.sitecoords(subtg,1),wdataset.sitecoords(subtg,2),10,'k','filled','MarkerFaceColor','g','Marker','+','MarkerEdgeColor','g'); hold on;
 
axis tight;
pdfwrite(['sitemap']);
