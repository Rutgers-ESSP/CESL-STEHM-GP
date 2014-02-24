function [f2s,sd2s,V2s,testlocs]=RegressHoloceneDataSets(dataset,testsitedef,modelspec,thetTGG,GISfp,ICE5G,trainsub,trainsubsubset,noiseMasks,testt,refyear)

%
% 
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 24 09:50:36 EST 2014



defval('testt',[-2000:10:2010]);
defval('refyear',2010);
defval('GISfp',[]);

if ~iscell(testt)
    testt=testt(:);
    testts=[];
else
    testts=testt;
end

cvfuncTGG=modelspec.cvfunc;
traincvTGG=modelspec.traincv;
if isfield(dataset,'bedrockMask')
    bedrockMask=dataset.bedrockMask;
else
    bedrockMask=@(a,b) [];
end

datid = dataset.datid;
istg = dataset.istg;
if isfield(dataset,'meantime')
    meantime=dataset.meantime;
else
    meantime=mean([dataset.time1(:) dataset.time2(:)],2);
end
lat=dataset.lat;
long=mod(dataset.long,360);
Y=dataset.Y;
Ycv=dataset.Ycv;
if isfield(dataset,'bedmsk')
    bedmsk=dataset.bedmsk;
else
    bedmsk=@(a,b) [];
end

limiting=dataset.limiting;
if isfield(dataset,'obsGISfp')
    obsGISfp=dataset.obsGISfp;
else
    obsGISfp=ones(size(lat));
end
compactcorr=dataset.compactcorr;
time1=dataset.time1;
time2=dataset.time2;
dY=dataset.dY;

testsites=testsitedef.sites;
testsites(:,3)=mod(testsites(:,3),360);
testnames=testsitedef.names;
testnames2=testsitedef.names2;

defval('trainsub',find((limiting==0).*(datid>0))); % only index points, and drop the global curve
defval('trainsubsubset',trainsub);

trainsub=intersect(trainsub,trainsubsubset);

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

dt1t1=dYears(meantime(trainsub),meantime(trainsub));
dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);

wcvfunc = @(x1,x2,thet) cvfuncTGG(x1,x2,dYears(x1,x2),thet,dy1y1,bedmsk(trainsub,trainsub),fp1fp1);

dt = abs(time2-time1)/4;

if max(dt)>0
    [dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),dY(trainsub),@(x1,x2) wcvfunc(x1,x2,thetTGG),2);
else
    dK=0; df=0; d2f=0; yoffset=0;
end

if iscell(GISfp)
    testsitefp = interp2(GISfp.long,GISfp.lat,GISfp.fp,testsites(:,3),testsites(:,2),'linear');
    testsitefp(find(testsites(:,2)>100))=1;
else
    testsitefp=ones(size(testsites,1));
end

testX=[];
testreg=[];
testfp=[];
for i=1:size(testsites,1)
    if iscell(testts)
        testt=testts{i}(:);
    end
	testX = [testX; ones(size(testt))*testsites(i,2) ones(size(testt))*testsites(i,3) testt];
	testreg = [testreg ; ones(size(testt))*testsites(i,1)];
	testfp=[testfp ; ones(size(testt))*testsitefp(i)];
end

Mref = eye(size(testX,1));
for i=1:size(testsites,1)
	
	sub1=find(testreg==testsites(i,1));
	sub2=intersect(sub1,find(testX(:,3)==refyear));
	
	Mref(sub1,sub2)=Mref(sub1,sub2)-1;

end

testGIAproju=zeros(size(testsites,1),1);
testGIAproj=zeros(size(testX,1),1);

for i=1:size(testsites,1)
	if testsites(i,1)>0
		testGIAproju(i)=interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsites(i,2),testsites(i,3));
		sub=find(testreg==testsites(i,1));
		testGIAproj(sub)=testGIAproju(i).*(testX(sub,3)-refyear);
	end
end

GIAproj=zeros(size(Y));

for jj=1:length(dataset.siteid)
    if dataset.siteid(jj)>0
        sub=find(datid==dataset.siteid(jj));       
        GIAproju(jj)=interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,lat(sub(jj)),long(sub(jj)));
        GIAproj(sub)=GIAproju(jj).*(meantime(sub)-refyear);
    end
end


X1 = [lat long meantime];

dt1t2 = dYears(X1(trainsub,3),testX(:,3));
dt2t2 = dYears(testX(:,3),testX(:,3));

dy1y2 = dDist(X1(trainsub,1:2),testX(:,1:2));
dy2y2 = dDist(testX(:,1:2),testX(:,1:2));

fp1fp2=bsxfun(@times,obsGISfp(trainsub)-1,testfp'-1)';

fp2fp2=bsxfun(@times,testfp-1,testfp'-1);


t1=X1(trainsub,3);
t2=testX(:,3);

testcv = @(thetas) cvfuncTGG(t1,t2,dt1t2,thetas,dy1y2,bedrockMask(datid(trainsub),testreg),fp1fp2);
testcv2 = @(thetas) cvfuncTGG(t2,t2,dt2t2,thetas,dy2y2,bedrockMask(testreg,testreg),fp2fp2) ;

dY2=sqrt(dY(trainsub).^2+dK);
Ycv2=Ycv(trainsub,trainsub)+diag(dK);

clear f2s V2s sd2s;
for i=1:size(noiseMasks,1)
	disp(i);
	[f2s(:,i),V2s(:,:,i),logp] = GaussianProcessRegression([],Y(trainsub)-GIAproj(trainsub)-yoffset,[],traincvTGG(t1,t1,dt1t1,thetTGG,Ycv2,dy1y1,bedmsk(trainsub,trainsub),fp1fp1),testcv(thetTGG.*noiseMasks(i,:))',testcv2(thetTGG.*noiseMasks(i,:)));
	sd2s(:,i)=sqrt(diag(V2s(:,:,i)));

end

parfor i=1:size(noiseMasks,1)
	if noiseMasks(i,8)==1
		f2s(:,i)=f2s(:,i)+testGIAproj;
	end
end

testlocs=testsitedef;
testlocs.X=testX;
testlocs.reg=testreg;
testlocs.fp=testfp;
testlocs.GIAproju=testGIAproju;
testlocs.GIAproj=testGIAproj;
