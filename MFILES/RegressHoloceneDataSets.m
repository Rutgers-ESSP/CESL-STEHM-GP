function [f2s,sd2s,V2s,testlocs,logp,passderivs,invcv]=RegressHoloceneDataSets(dataset,testsitedef,modelspec,thetTGG,trainsub,noiseMasks,testt,refyear,collinear,passderivs,invcv)

% [f2s,sd2s,V2s,testlocs,logp,passderivs,invcv]=RegressHoloceneDataSets(dataset,
%       testsitedef,modelspec,thetTGG,[trainsub],noiseMasks,testt,[refyear],
%       [collinear],[passderivs],[invcv])
%
% INPUT:
%    dataset - sea-level data structure dataset (with fields datid, istg,
%              time1, time2, lat, long, Y, Ycv, limiting, compactcorr, and dY)
%    testsitedef - site definition structure (with fields sites, names, names2,
%                  and optionally firstage and GIA, where sites is a N x 3
%                   matrix with columns corresponding to siteid, lat, and long)
%    modelspec - model specification structure (with fields cvfunc and traincv,
%                and optionally derivatives dcvfunc and ddcvfunc)
%    thetTGG - hyperparameters
%    trainsub - training subset indices
%    noiseMasks - matrix of noise masks, which will be multiplied row-wise with
%                 thetTGG
%    testt - target ages
%    refyear - reference year for GIA calculations (default = 2010)
%    collinear - column of thetTGG corresponding to the amplitude of the linear
%                rate term; if equal to zero or noiseMasks(i,collinear) = 1, the
%                GIA rate will be subtracted before regression and then added
%                back in; otherwise it will not be added back in
%    passderivs - derivatives passed along to speed calculation
%    invcv - inverse of covariance matrix, passed along to speed calculations
%
% OUTPUT:
%    f2s - estimated levels (columns corresponding to rows of noiseMasks)
%    sd2s - estimated standard deviaitons
%    V2s - estimated covariance
%    testlocs - site definition structure
%    logp - log likelihood
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Oct 26 13:43:00 EDT 2014


defval('testt',[-2010:20:2010]);
defval('refyear',2010);
defval('collinear',0);
defval('passderivs',[]);
defval('invcv',[]);

if ~iscell(testt)
    testt=testt(:);
    testts=[];
else
    testts=testt;
end

cvfuncTGG=modelspec.cvfunc;
traincvTGG=modelspec.traincv;

analyticalderiv = 0;
if isfield(modelspec,'ddcvfunc')
    if length(modelspec.ddcvfunc)>0
        dcvfuncTGG = modelspec.dcvfunc;
        ddcvfuncTGG = modelspec.ddcvfunc;
        analyticalderiv = 1;
    end
end


datid = dataset.datid;
istg = dataset.istg;
if isfield(dataset,'meantime')
    meantime=dataset.meantime;
else
    meantime=mean([dataset.time1(:) dataset.time2(:)],2);
end
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y;
Ycv=dataset.Ycv;

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

% compaction correction
dY=sqrt(dY.^2+(compactcorr.*thetTGG(end)).^2);
Ycv=Ycv+diag(compactcorr.*thetTGG(end)).^2;

testsites=testsitedef.sites;
testnames=testsitedef.names;
testnames2=testsitedef.names2;
if isfield(testsitedef,'firstage')
    firstage=testsitedef.firstage;
else
    firstage=-Inf*ones(size(testsites));
end

defval('trainsub',find((limiting==0).*(datid>0))); % only index points, and drop the global curve

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

dt1t1=dYears(meantime(trainsub),meantime(trainsub));
dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);

dK=0; df=0; d2f=0; yoffset=0;
dt = abs(time2-time1)/4;
if length(passderivs)>0
    dK = passderivs.dK;
    yoffset = passderivs.yoffset;
else
    
    if max(dt)>0
        if analyticalderiv
            mspec.cvfunc = @(x1,x2,xxx,yyy) cvfuncTGG(x1,x2,dYears(x1,x2),thetTGG,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
            mspec.dcvfunc = @(x1,x2,xxx,yyy) dcvfuncTGG(x1,x2,dYears(x1,x2),thetTGG,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
            mspec.ddcvfunc = @(x1,x2,xxx,yyy) ddcvfuncTGG(x1,x2,dYears(x1,x2),thetTGG,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
            [dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),dY(trainsub),mspec,1,[1:length(trainsub)]');
        else
            wcvfunc = @(x1,x2,thet,xxx,yyy) cvfuncTGG(x1,x2,dYears(x1,x2),thet,dy1y1(xxx,yyy)',fp1fp1(xxx,yyy)');
           
           [dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),dY(trainsub),@(x1,x2,r1,r2) wcvfunc(x1,x2,thetTGG,r1,r2),1,[1:length(trainsub)]');
        end
        
    end
    passderivs.dK = dK;
    passderivs.yoffset = yoffset;
    passderivs.df = df;
    passderivs.d2d = d2f;
end

if isfield(testsitedef,'GISfp')
    testsitefp=testsitedef.GISfp;
else
    testsitefp=ones(size(testsites(:,1)));
end

if isfield(testsitedef,'GIA')
    testGIAproju=testsitedef.GIA;
else
    testGIAproju=zeros(size(testsites(:,1)));
end

testX=[];
testreg=[];
testfp=[];
for i=1:size(testsites,1)
    if iscell(testts)
        testt=testts{i}(:);
    end
    subtimes=find(testt>=firstage(i));
    testX = [testX; ones(size(testt(subtimes)))*testsites(i,2) ones(size(testt(subtimes)))*testsites(i,3) testt(subtimes)];
    testreg = [testreg ; ones(size(testt(subtimes)))*testsites(i,1)];
    testfp=[testfp ; ones(size(testt(subtimes)))*testsitefp(i)];
end

testGIAproj=zeros(size(testX,1),1);

for i=1:size(testsites,1)
    sub=find(testreg==testsites(i,1));
    testGIAproj(sub)=testGIAproju(i).*(testX(sub,3)-refyear);
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

testcv = @(thetas) cvfuncTGG(t1,t2,dt1t2,thetas,dy1y2,fp1fp2);
testcv2 = @(thetas) cvfuncTGG(t2,t2,dt2t2,thetas,dy2y2,fp2fp2) ;

dY2=sqrt(dY(trainsub).^2+dK);
Ycv2=Ycv(trainsub,trainsub)+diag(dK);

clear f2s V2s sd2s;
for i=1:size(noiseMasks,1)
    disp(i);
    [f2s(:,i),V2s(:,:,i),logp,~,~,~,invcv] = GaussianProcessRegression([],Y(trainsub)-yoffset,[],traincvTGG(t1,t1,dt1t1,thetTGG,Ycv2,dy1y1,fp1fp1),testcv(thetTGG.*noiseMasks(i,:))',testcv2(thetTGG.*noiseMasks(i,:)),invcv);
    sd2s(:,i)=sqrt(diag(V2s(:,:,i)));
    if (collinear==0)||(collinear>size(noiseMasks,2))
        f2s(:,i)=f2s(:,i)+testGIAproj;
    elseif noiseMasks(i,collinear)==1
        f2s(:,i)=f2s(:,i)+testGIAproj;
    end
    
end

testlocs=testsitedef;
testlocs.X=testX;
testlocs.reg=testreg;
testlocs.fp=testfp;
testlocs.GIAproju=testGIAproju;
testlocs.GIAproj=testGIAproj;

