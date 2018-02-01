% Identify sites for prediction
% 
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-01-30 23:39:47 -0500

clear oldest oldestnear youngest testsitedef;
wdataset = datasets{1};
distfrom = dDist(wdataset.sitecoords,wdataset.sitecoords);
for ii=1:length(wdataset.sitenames)
    sub = find(wdataset.datid==wdataset.siteid(ii));
    if length(sub)>0
        oldest(ii)=min(union(wdataset.time1(sub),wdataset.time2(sub)));
        youngest(ii)=max(union(wdataset.time1(sub),wdataset.time2(sub)));
    else
        oldest(ii)=2000;
        youngest(ii)=2014;
    end
end
for ii=1:length(wdataset.sitenames)
    sub=find(distfrom(ii,:)<.5);
    if wdataset.sitecoords(ii,1)>360
        sub=ii;
    end
    oldestnear(ii) = min(oldest(sub));
end

clear testsitedef;
testsitedef.sites=[0 1e6 1e6];
testsitedef.names={'GSL'};
testsitedef.names2={'GSL'};
testsitedef.firstage=min(oldest);
testsitedef.oldest=min(oldest);
testsitedef.youngest=2014;

testsitedefGSL=testsitedef;

dosub=find(wdataset.siteid>0);
avail=ones(length(wdataset.siteid),1); % make more restrictive if you want to exclude some data sets
wsitelen=wdataset.sitelen;
wsitelen(find(wdataset.siteid<9000))=wsitelen(find(wdataset.siteid<9000))/10;
[s,si]=sort(wsitelen(dosub));
dosub=dosub(si(end:-1:1));
dosub=dosub([find(wdataset.siteid(dosub)>9000)' find(wdataset.siteid(dosub)<=9000)']); 

for iii=1:length(dosub)
    ii=dosub(iii);
    if avail(ii)
        sub=find((wdataset.datid==wdataset.siteid(ii)));
        if length(sub)>0
            testsitedef.sites(end+1,:)=mean([wdataset.datid(sub) wdataset.lat(sub) wdataset.long(sub)],1);
            wdist = dDist(testsitedef.sites(end,2:3),wdataset.sitecoords);
            avail=avail.*(wdist>.2); % keep sites at least .2 degrees apart
            testsitedef.names2={testsitedef.names2{:}, wdataset.sitenames{ii}};
            
            sublett=setdiff(1:length(wdataset.sitenames{ii}),strfind(wdataset.sitenames{ii},' '));
            testsitedef.names={testsitedef.names{:}, wdataset.sitenames{ii}(sublett)};
            testsitedef.firstage = [testsitedef.firstage oldestnear(ii)];
            testsitedef.oldest = [testsitedef.oldest oldest(ii)];
            testsitedef.youngest = [testsitedef.youngest youngest(ii)];
        end
    end
    
end

testsitedef.GIA=zeros(size(testsitedef.sites(:,2)));
if exist('ICE5G','var')
    ICE5G.lat=ICE5Glat;
    ICE5G.long=ICE5Glon;
    ICE5G.gia=ICE5Ggia;

    if length(ICE5G.gia)>1
        testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
        testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;
        testsitedef.GIA(find(isnan(testsitedef.GIA)))=0;
    end
    else
    ICE5G=0;
end