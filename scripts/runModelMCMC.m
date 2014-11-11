% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Nov 11 00:16:55 EST 2014

% now set up and run parallel tempering MCMC

ii=1;
subfixed=modelspec(trainspecs(ii)).subfixed;
trainsub=trainsubsubset{ii};
wdataset=datasets{trainsets(ii)};

subnotfixed=setdiff(1:length(thetTGG{ii})-1,subfixed);
Mfixed=sparse(length(subnotfixed),length(thetTGG{ii})-1);
for i=1:length(subnotfixed)
Mfixed(i,subnotfixed(i))=1;
end
fixedvect = 0*thetTGG{ii}(1:end-1); fixedvect(subfixed)=thetTGG{ii}(subfixed);

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

            dt1t1=dYears(wdataset.meantime(trainsub),wdataset.meantime(trainsub));
            dy1y1 = dDist([wdataset.lat(trainsub) wdataset.long(trainsub)],[wdataset.lat(trainsub) wdataset.long(trainsub)]);
            fp1fp1=bsxfun(@times,wdataset.obsGISfp(trainsub)-1,wdataset.obsGISfp(trainsub)'-1);


mspec.cvfunc = @(x1,x2,thet,xxx,yyy) modelspec(trainspecs(ii)).cvfunc(x1,x2,dYears(x1,x2),thet(1:end-1)*Mfixed+fixedvect,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
mspec.dcvfunc = @(x1,x2,thet,xxx,yyy) modelspec(trainspecs(ii)).dcvfunc(x1,x2,dYears(x1,x2),thet(1:end-1)*Mfixed+fixedvect,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
mspec.ddcvfunc = @(x1,x2,thet,xxx,yyy) modelspec(trainspecs(ii)).ddcvfunc(x1,x2,dYears(x1,x2),thet(1:end-1)*Mfixed+fixedvect,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');

maxcompactfactor=1;
dolb = [ modelspec(ii).lb(subnotfixed) 1e-6];
doub = [modelspec(ii).ub(subnotfixed) maxcompactfactor];
subnotfixed=union(subnotfixed,length(thetTGG{ii}));

[thet,logp,acceptedcount]=SLNIGPSamplePT(wdataset.meantime(trainsub),wdataset.Y(trainsub), ...
                                       wdataset.dt(trainsub),wdataset.dY(trainsub),@(thet) ...
                                       diag(thet(end)*wdataset.compactcorr(trainsub)).^2, mspec, ...
                                       thetTGG{ii}(subnotfixed),dolb,doub,[1:length(trainsub)]',8,10000,1,1,.1,'samplesPT','samplesPT');
