function [TGdata2,TGdata,thetL,TGmodellocal] = GPSmoothTideGauges(targcoords,addlsites,thetL0,winlength,thinlength,minlength,optimizemode,psmsldir,gslfile,ROOTDIR)

% Find tide gauge sites to include, based on length and proximity criteria, then
% fit GP model to them in order to interpolate and take running averge.
%
% Should tweak to use proxy data set rather than prescribed list of sites to look for 
% tide gauges near.
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Feb 28 10:23:00 EST 2014
%

defval('minlength',[120 70 20]);  % minimum length for global curves, most purposes and minimum
defval('thetL0',[]);
defval('winlength',11);
defval('thinlength',winlength-1);
defval('GIAanchoryear',2005);
defval('optimizemode',1.1); % set to 1.0 for local optimization only

defval('targcoords',[
34.97  -76.38; %Tump Point, NC
35.89  -75.64; %Sand Point, NC
39.085 -74.811; % Cape May Courthouse, NJ
39.50  -74.415; % Leeds Point, NJ
30.6   -81.7; % Nassau, FL
44.72  -63.23; % Chezzetcook, Nova Scotia
]);

defval('addlsites',[]);

defval('ROOTDIR','~/Dropbox/Consulting/Risky Business/Code/RiskyBusinessScience/slr');
addpath(fullfile(ROOTDIR,'MFILES'));
IFILES=fullfile(ROOTDIR,'../../IFILES/slr/');

%defval('giafile',fullfile(IFILES,'dsea250.1grid.ICE5Gv1.3_VM2_L90_2012.nc'));
defval('psmsldir',fullfile(IFILES,'rlr_annual'));
defval('gslfile',fullfile(IFILES,'CSIRO_Recons_gmsl_yr_2011.csv'));

%giamodel.gia=ncread(giafile,'Dsea_250');
%giamodel.lat=ncread(giafile,'Lat');
%giamodel.long=ncread(giafile,'Lon');

[TGcoords,TGrsl,TGrslunc,TGid,TGsiteid,sitenames,TGsitecoords,sitelen]=ReadPSMSLData(1,1000,minlength(3),psmsldir,gslfile);
sub1=find((sitelen>minlength(1)));
sub1=union(sub1,find(TGsiteid==0));

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

sub2=[];
for ii=1:size(targcoords,1)
    TGdist=dDist(TGsitecoords,targcoords(ii,:));
    [m,mi]=min(TGdist);
    mi=union(mi,find((TGdist<2).*(sitelen'>minlength(2))));
    sub2=union(sub2,mi);
end

addlsites=addlsites(:)';
subaddlsites=find(ismember(TGsiteid,addlsites));

sitesub=union(sub1,sub2);
sitesub=union(sitesub,subaddlsites);
sub=find(ismember(TGid,TGsiteid(sitesub)));

clear TGdata;
TGdata.datid=TGid(sub);
TGdata.time1=TGcoords(sub,3);
TGdata.time2=TGcoords(sub,3);
TGdata.meantime=TGcoords(sub,3);
TGdata.limiting=zeros(size(TGid(sub)));
TGdata.Y=TGrsl(sub);
TGdata.dY =TGrslunc(sub);
TGdata.compactcorr=zeros(size(TGid(sub)));
TGdata.istg = ones(size(TGid(sub)));
TGdata.lat=TGcoords(sub,1);
TGdata.long=TGcoords(sub,2);;
TGdata.Ycv=sparse(diag(TGrslunc(sub).^2));
TGdata.siteid=TGsiteid(sitesub);
TGdata.sitenames=sitenames(sitesub);
TGdata.sitecoords=TGsitecoords(sitesub,:);
TGdata.sitelen=sitelen(sitesub);

%%%%%%%%%%

GPSLDefineCovFuncs;

TGmodellocal.cvfunc=@(t1,t2,dt1t2,thetas,dy1y2,fp1fp2) kDP(t1,t2,thetas(1)) + kMatG(dt1t2,thetas(2:4)) + kDELTAG(dy1y2,thetas(5));
TGmodellocal.traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) TGmodellocal.cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;


tluL = [
1   0.1    100 % DP

5   0.1   1e3   % MatG
1   .1  300
1.5  .5  5.5

1  0  1e4 % offset

];

TGmodellocal.thet0=tluL(:,1)';
TGmodellocal.lb=tluL(:,2)';
TGmodellocal.ub=tluL(:,3)';
TGmodellocal.subfixed=[];
TGmodellocal.sublength=[];

clear TGdata2 thetL;
TGdata2.datid=[];
TGdata2.time1=[];
TGdata2.Y=[];
TGdata2.dY =[];
TGdata2.lat=[];
TGdata2.long=[];
TGdata2.Ycv=[];

for nn=1:length(TGdata.siteid)

    targsite=TGdata.siteid(nn);
    sub=find(TGdata.datid==targsite);
    subS=find(TGdata.siteid==targsite);
    disp(TGdata.sitenames{subS});
    
    TGdatasub=TGdata;
    shortenfields={'datid','time1','time2','meantime','limiting','Y','dY','compactcorr','istg','lat','long'};
    for jj=1:length(shortenfields)
        TGdatasub.(shortenfields{jj}) =  TGdatasub.(shortenfields{jj})(sub);
    end
    shortenfields={'siteid','sitenames','sitelen'};
    for jj=1:length(shortenfields)
        TGdatasub.(shortenfields{jj}) =  TGdatasub.(shortenfields{jj})(subS);
    end
    TGdatasub.sitecoords=TGdatasub.sitecoords(subS,:);
    TGdatasub.Ycv=sparse(diag(TGdatasub.dY.^2));
    TGdatasub.Ycv0=sparse(diag(TGdatasub.dY.^2));
    if length(thetL0)==0
        [thetL(nn,:)]=OptimizeHoloceneCovariance(TGdatasub,TGmodellocal,optimizemode)
    elseif size(thetL0,1)>1
        thetL(nn,:)=thetL0(nn,:);
    else
        thetL(nn,:)=thetL0;
    end

    clear TGtestsitedef;
    TGtestsitedef.sites=[TGdatasub.siteid TGdatasub.sitecoords];
    TGtestsitedef.names=TGdatasub.sitenames;
    TGtestsitedef.names2=TGtestsitedef.names;
    TGtestsitedef.t = floor(min(TGdata.time1(sub)-winlength/2)):ceil(max(TGdata.time1(sub)+winlength/2));
    noiseMasks = ones(1,size(thetL,2));
    noiseMasklabels={'full'};
    trainsub=find(TGdatasub.limiting==0);
    [TGf,TGsd,TGV,TGtestlocs]=RegressHoloceneDataSets(TGdatasub,TGtestsitedef,TGmodellocal,thetL(nn,:),[],[],trainsub,[],noiseMasks,TGtestsitedef.t,GIAanchoryear);
    
    % check for bad fit, and do a global optimization if need me
    if min(diag(TGV))<0
        if (length(thetL0)==0) .* (optimizemode<1.1)
             [thetL(nn,:)]=OptimizeHoloceneCovariance(TGdatasub,TGmodellocal,1.1)
             [TGf,TGsd,TGV,TGtestlocs]=RegressHoloceneDataSets(TGdatasub,TGtestsitedef,TGmodellocal,thetL(nn,:),[],[],trainsub,[],noiseMasks,TGtestsitedef.t,GIAanchoryear);
        end
    end

    % now take running averages
    Mop=abs(bsxfun(@minus,TGtestlocs.X(:,3),TGtestlocs.X(:,3)'))<(winlength/2);
    Mop=Mop.*bsxfun(@eq,TGtestlocs.reg,TGtestlocs.reg');
    adder=sum(Mop,2);
    sub=find(adder==winlength);
    Mop=Mop(sub,:)/winlength;
    t2=Mop*TGtestlocs.X(:,3);
    reg2=round(Mop*TGtestlocs.reg*10)/10;
    selsub=[];
    for ii=1:length(TGdatasub.siteid)
        subO=find(TGdatasub.datid==TGdatasub.siteid(ii));
        sub=find((reg2==TGdatasub.siteid(ii)).*(t2>=min(TGdatasub.time1(subO))).*(t2<=max(TGdatasub.time1(subO))));
        sub3=sub(1:thinlength:end);
        selsub=[selsub ; sub3];
    end
    Mop=Mop(selsub,:);

    TGf2=Mop*TGf;
    TGV2=Mop*TGV*Mop';
    TGsd2=sqrt(diag(TGV2));

    TGdata2.datid=[TGdata2.datid ; Mop*TGtestlocs.reg];
    TGdata2.time1=[TGdata2.time1 ; Mop*TGtestlocs.X(:,3)];
    TGdata2.Y=[TGdata2.Y ; TGf2];
    TGdata2.dY =[TGdata2.dY ; TGsd2];
    TGdata2.lat=[TGdata2.lat ; Mop*TGtestlocs.X(:,1)];
    TGdata2.long=[TGdata2.long ; Mop*TGtestlocs.X(:,2)];
    TGdata2.Ycv(length(TGdata2.Ycv)+[1:length(TGf2)],length(TGdata2.Ycv)+[1:length(TGf2)])=sparse(TGV2);

    TGdata2.compactcorr=zeros(size(TGdata2.datid));
    TGdata2.time2=TGdata2.time1;
    TGdata2.meantime=TGdata2.time1;
    TGdata2.limiting=zeros(size(TGdata2.datid));
    TGdata2.istg = ones(size(TGdata2.datid));
    TGdata2.siteid=TGdata.siteid;
    TGdata2.sitenames=TGdata.sitenames;
    TGdata2.sitecoords=TGdata.sitecoords;
    TGdata2.sitelen=TGdata.sitelen;

end
