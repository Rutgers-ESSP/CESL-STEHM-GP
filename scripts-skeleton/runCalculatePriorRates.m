% Calculate rates from priors
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Jan 04 16:50:49 EST 2016


%% Determine prior rates

% define intervals focused upon
firstyearspriors=[0    0  300 700 1000 1400 1400 1600 1800 1860 1900];
lastyearspriors= [700 300 700 1000 1400 1800 1600 1800 1900 1900 2000];
qlevsmmpriors=[.05 .167 .5 .833 .95];


clear priorslopef priorslopesd priorslopedf priormxtomnq;
for qqq=1:length(thetTGG)
    ms = modelspec(trainspecs(qqq));
    wdat = datasets{trainsets(qqq)};
    dothet=thetTGG{qqq};
    
    % create noiseMasks -- set to zero only for amplitudes of
    % terms you want to exclude; if you set to zero for a non-amplitude
    % term it is likely to create problems
    
    noiseMasks = ones(1,length(thetTGG{qqq}));
    noiseMasks(1,[ ms.subampnoise ]  )=0; %without noise
    noiseMasklabels={'denoised'};
    collinear=ms.subamplinear(1);


    flattenersub=find((wdat.datid==0).*(wdat.dY>1e3));

    if length(flattenersub)==0
        datsub=find(wdat.datid==0); datsub=datsub(end);
        sitesub=find(wdat.siteid==0);
        wd=SubsetDataStructure(wdat,datsub,sitesub);
        wd.dY = 100000;
        wd.Ycv = wd.dY^2;
        wd.meantime=100000;
        wd.time1=wd.meantime; wd.time2=wd.meantime;
    else
          datsub=flattenersub;
          sitesub=find(wdat.siteid==0);
          wd=SubsetDataStructure(wdat,datsub,sitesub);

    end

    clear wdef;
    wdef.sites=[0 1e6 1e6];
    wdef.names={'GSL'};
    wdef.names2={'GSL'};
    wdef.firstage=0;
    wdef.oldest=0;
    wdef.youngest=2014;
    trainsub=find(wd.limiting==0);
    wtestt=testt;
    [wf,wsd,wV,wloc]=RegressHoloceneDataSets(wd,wdef,ms,dothet,trainsub,noiseMasks(1,:),wtestt,refyear,collinear);        
    [wfslope,wsdslope,wfslopediff,wsdslopediff,wdiffplus,wdiffless]=SLRateCompare(wf,wV,wloc.sites,wloc.reg,wloc.X(:,3),firstyearspriors,lastyearspriors);
    priorslopef(qqq,:)=wfslope(1,:);
    priorslopesd(qqq,:)=wsdslope(1,:);
    
    samps=lhsnorm(wf(:,1),wV(:,:,1),1000);
    sub2=find((wloc.X(:,3)<1900).*(wloc.X(:,3)>=0));
    mxtomn=max(samps(:,sub2),[],2)-min(samps(:,sub2),[],2);
    priormxtomnq(qqq,:) = quantile(mxtomn,qlevsmmpriors);
end
