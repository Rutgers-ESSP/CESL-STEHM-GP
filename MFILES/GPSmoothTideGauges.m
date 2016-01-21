function [TGdata2s,thetL,TGmodellocal] = GPSmoothTideGauges(TGdata,winlengths,optimizemode,thinlengths,thetL0,thinyrstart,noiseMask)

% TGdata2 = GPSmoothTideGauges(TGdata,[winlength],[optimizemode],[thinlength],[thetL0])
%
% Fit GP model to TGdata site-by-site in order to interpolate and take running averge.
%
%   TGdata: data structure with tide gauge data to smooth
%   winlength: window length for running averages, or list of window lengths (default = 11)
%   optimizemode: optimization mode (default = 1.1 for global optimization, set to 1.0 for local)
%   thinlength: spacing of data to output, same length as winlength (default = winlength - 1)
%   thetL0: hyperparameters of GP model (default optimizes)
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Dec 24 11:53:13 EST 2015
%

defval('thetL0',[]);
defval('winlengths',11);
defval('thinlengths',max(1,winlengths));
defval('optimizemode',1.1); % set to 1.0 for local optimization only
defval('thinyrstart',[]);
defval('noiseMask',[]);

CESLDefineCovFuncs;

TGmodellocal.cvfunc=@(t1,t2,dt1t2,thetas,dy1y2,fp1fp2) kDP(t1,t2,thetas(1)) + kMatG(dt1t2,thetas(2:4)) + kDELTA(dt1t2,thetas(5)) + kCONST(thetas(6));
TGmodellocal.traincv = @(t1,t2,dt1t2,thetas,errcv,ad,fp1fp2) TGmodellocal.cvfunc(t1,t2,dt1t2,thetas,ad,fp1fp2) + errcv;

tluL = [
1   0.1    100 % DP

5   0.1   1e3   % MatG
10   1  300
1.5  .5  5.5

    
5  0.1 1e3 % annual variability    
1  0  1e4 % offset

];

TGmodellocal.thet0=tluL(:,1)';
TGmodellocal.lb=tluL(:,2)';
TGmodellocal.ub=tluL(:,3)';
TGmodellocal.subfixed=[];
TGmodellocal.sublength=[];

clear TGdata2s thetL;

for kkk=1:length(winlengths)
    TGdata2s{kkk}.datid=[];
    TGdata2s{kkk}.time1=[];
    TGdata2s{kkk}.Y=[];
    TGdata2s{kkk}.dY =[];
    TGdata2s{kkk}.lat=[];
    TGdata2s{kkk}.long=[];
    TGdata2s{kkk}.Ycv=[];
end
        
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
    if length(thinyrstart)==0
        TGtestsitedef.t = floor(min(TGdata.time1(sub)-max(winlengths)/2)):min([thinlengths 1]):ceil(max(TGdata.time1(sub)+max(winlengths)/2));
    else
        TGtestsitedef.t = floor(min(min(TGdata.time1(sub)-max(winlengths)/2),thinyrstart)):min([thinlengths 1]):ceil(max(TGdata.time1(sub)+max(winlengths)/2));
    end
    
    if length(noiseMask)==0
        noiseMasks = ones(1,size(thetL,2));
    else
        noiseMasks=noiseMask;
    end
    
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

    for kkk=1:length(winlengths)
    
        winlength=winlengths(kkk);
        thinlength=thinlengths(kkk);

        % now take running averages
        Mop=abs(bsxfun(@minus,TGtestlocs.X(:,3),TGtestlocs.X(:,3)'))<=(winlength/2);
        Mop=Mop.*bsxfun(@eq,TGtestlocs.reg,TGtestlocs.reg');
        adder=sum(Mop,2);
        sub=find(adder>=winlength);
        Mop=bsxfun(@rdivide,Mop(sub,:),adder(sub));
        t2=Mop*TGtestlocs.X(:,3);
        t2=round(t2*1000)/1000;
        if length(thinyrstart)==0
            [selt2, selsub] = intersect(t2,min(TGdatasub.time1):thinlength:max(TGdatasub.time1));
        else
            timegrid = thinyrstart:thinlength:(max(TGdatasub.time1)+thinlength);
            subq = find((timegrid>=min(TGdatasub.time1)).*(timegrid<=max(TGdatasub.time1)));
            if subq(1)>1
                if length(find(timegrid==min(TGdatasub.time1)))==0
                    subq = [subq(1)-1 subq];
                end
            end
            if subq(end)<length(timegrid)
                if length(find(timegrid==max(TGdatasub.time1)))==0
                    subq = [subq subq(end)+1];
                end
            end
            
            [selt2, selsub] = intersect(t2,timegrid(subq));
        end
        
        Mop=Mop(selsub,:);
        t2=t2(selsub);

        TGf2=Mop*TGf;
        TGV2=Mop*TGV*Mop';
        TGsd2=sqrt(diag(TGV2));

        % and put back in a data structure
        TGdata2s{kkk}.datid=[TGdata2s{kkk}.datid ; round(Mop*TGtestlocs.reg)];
        TGdata2s{kkk}.time1=[TGdata2s{kkk}.time1 ; t2];
        TGdata2s{kkk}.Y=[TGdata2s{kkk}.Y ; TGf2];
        TGdata2s{kkk}.dY =[TGdata2s{kkk}.dY ; TGsd2];
        TGdata2s{kkk}.lat=[TGdata2s{kkk}.lat ; Mop*TGtestlocs.X(:,1)];
        TGdata2s{kkk}.long=[TGdata2s{kkk}.long ; Mop*TGtestlocs.X(:,2)];
        TGdata2s{kkk}.Ycv(length(TGdata2s{kkk}.Ycv)+[1:length(TGf2)],length(TGdata2s{kkk}.Ycv)+[1:length(TGf2)])=sparse(TGV2);

        TGdata2s{kkk}.compactcorr=zeros(size(TGdata2s{kkk}.datid));
        TGdata2s{kkk}.time2=TGdata2s{kkk}.time1;
        TGdata2s{kkk}.meantime=TGdata2s{kkk}.time1;
        TGdata2s{kkk}.limiting=zeros(size(TGdata2s{kkk}.datid));
        TGdata2s{kkk}.istg = ones(size(TGdata2s{kkk}.datid));
        TGdata2s{kkk}.siteid=TGdata.siteid;
        TGdata2s{kkk}.sitenames=TGdata.sitenames;
        TGdata2s{kkk}.sitecoords=TGdata.sitecoords;
        TGdata2s{kkk}.sitelen=TGdata.sitelen;
    end
    



end
    
if length(winlengths)==1
    TGdata2s = TGdata2s{1};
end
