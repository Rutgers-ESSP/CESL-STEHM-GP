% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 30 17:07:28 EST 2014

% extract the detrended North Carolina sea-level curve

sitesub=find(strcmpi('North Carolina-Sand Point',wdataset.sitenames));
sitesub=union(sitesub,find(strcmpi('North Carolina-Tump Point',wdataset.sitenames)));
datsub=find(ismember(wdataset.datid,wdataset.siteid(sitesub)));

nmms=ones(size(thetNC)); nmms(3)=0;

testtNC = [-1000:20:2000 2010];

wtestsitedef=testsitedefGSL;
wtestsitedef.sites(:,2:3)=mean(wdataset.sitecoords(sitesub,:),1);
[NC_pseudoGSL,NC_pseudoGSLsd,NC_pseudoGSLV,NClocs]=RegressHoloceneDataSets(wdataset,testsitedefGSL,ms,thetNC,datsub,nmms,testtNC,refyear,3);
NC_yrs=NClocs.X(:,3);

[NC_pseudoGSL,NC_pseudoGSLV,NC_pseudoGSLsd]=DetrendSLReconstruction(NC_pseudoGSL,NC_pseudoGSLV+ eye(size(NC_pseudoGSLV))*50^2,NClocs.sites,NClocs.reg,NC_yrs,0,1800,2000);
Mrefms=eye(length(NC_pseudoGSL));
subms=find((NC_yrs<=1800).*(NC_yrs>=0));
Mrefms(:,subms)=Mrefms(:,subms)-1/length(subms);
NC_pseudoGSL=Mrefms*NC_pseudoGSL;
NC_pesudoGSLsd=sqrt(diag(Mrefms*NC_pseudoGSLV*Mrefms'));
