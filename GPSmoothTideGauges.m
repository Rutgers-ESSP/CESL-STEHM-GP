function TGdata2 = GPSmoothTideGauges(TGdata,winlength,optimizemode,thinlength,thetL0)

% TGdata2 = GPSmoothTideGauges(TGdata,[winlength],[optimizemode],[thinlength],[thetL0])
%
% Fit GP model to TGdata site-by-site in order to interpolate and take running averge.
%
%   TGdata: data structure with tide gauge data to smooth
%   winlength: window length for running averages (default = 11)
%   optimizemode: optimization mode (default = 1.1 for global optimization, set to 1.0 for local)
%   thinlength: spacing of data to output (default = winlength - 1)
%   thetL0: hyperparameters of GP model (default optimizes)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue May 20 17:59:24 EDT 2014
%

defval('thetL0',[]);
defval('winlength',11);
defval('thinlength',winlength-1);
defval('optimizemode',1.1); % set to 1.0 for local optimization only

TGDefineCovFuncs;

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

    % find the subset of data for each site
    targsite=TGdata.siteid(nn);
    sub=find(TGdata.datid==targsite);
    subS=find(TGdata.siteid==targsite);
    disp(TGdata.sitenames{subS});    
    TGdatasub=SubsetDataStructure(TGdata,sub,subS);

    % optimize GP hyperparameters unless specified
    if length(thetL0)==0
        [thetL(nn,:)]=OptimizeHoloceneCovariance(TGdatasub,TGmodellocal,optimizemode);
    elseif size(thetL0,1)>1
        thetL(nn,:)=thetL0(nn,:);
    else
        thetL(nn,:)=thetL0;
    end

    % predict tide gauge levels
    clear TGtestsitedef;
    TGtestsitedef.sites=[TGdatasub.siteid TGdatasub.sitecoords];
    TGtestsitedef.names=TGdatasub.sitenames;
    TGtestsitedef.names2=TGtestsitedef.names;
    TGtestsitedef.t = floor(min(TGdata.time1(sub)-winlength/2)):ceil(max(TGdata.time1(sub)+winlength/2));
    noiseMasks = ones(1,size(thetL,2));
    noiseMasklabels={'full'};
    trainsub=find(TGdatasub.limiting==0);
    [TGf,TGsd,TGV,TGtestlocs]=RegressHoloceneDataSets(TGdatasub,TGtestsitedef,TGmodellocal,thetL(nn,:),trainsub,noiseMasks,TGtestsitedef.t);
    
    % check for bad fit, and do a global optimization if need me
    if min(diag(TGV))<0
        if (length(thetL0)==0) .* (optimizemode<1.1)
             [thetL(nn,:)]=OptimizeHoloceneCovariance(TGdatasub,TGmodellocal,1.1)
             [TGf,TGsd,TGV,TGtestlocs]=RegressHoloceneDataSets(TGdatasub,TGtestsitedef,TGmodellocal,thetL(nn,:),trainsub,noiseMasks,TGtestsitedef.t);
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

    % and put back in a data structure
    TGdata2.datid=[TGdata2.datid ; round(Mop*TGtestlocs.reg)];
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
