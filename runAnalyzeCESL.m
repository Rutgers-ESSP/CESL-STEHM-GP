% Master script for Common Era proxy + tide gauge sea-level analysis
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Nov 04 18:09:17 EST 2014

dosldecomp = 0;

pd=pwd;
%addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pd '/MFILES']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='141103';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);

GIAfiles=([pd '/../GIA/RSL4/rsl*.out.gz']);

%
firsttime=-1000;

runSetupHoloceneCovariance;
runImportHoloceneDataSets;

save(savefile,'datasets','modelspec');

% let's start with full data set, but no sites requiring compaction correction 
% then add a compaction correction as a last step

%trainsets=[3 4 1 5];
%trainspecs=[1 1 1 1];
%trainlabels={'trTGPXnoEH','trTGPXGSLnoEH','trTGPX','trTGPXGSL'};

% need to see what happens if we drop EH from a constrained data set 


trainsets=[]; trainspecs=[]; trainlabels={};
             

% $$$ for ii=1:length(modelspec)
% $$$     trainsets(ii) = 5;
% $$$     trainspecs(ii) = ii;
% $$$     trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(ii).label]};
% $$$ end

% we will run: GLMW (5), LMW (2), GLW (4), GMW (9), LW (7), GW (8), MW (10)

trainspecs = [5 2 4 9 7 8 10];
trainsets = [5*ones(size(trainspecs)) 1*ones(size(trainspecs))];
trainspecs = [trainspecs trainspecs];

for ii=1:length(trainsets)
    trainlabels = {trainlabels{:}, [datasets{trainsets(ii)}.label '_' modelspec(trainspecs(ii)).label]};
end

trainrange=[100 100 2000 2000];
trainfirsttime = -1000;

clear thetTGG thethist trainsubsubset logp;
for ii=1:length(trainspecs)
    % first only fit ones without a compaction correction
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thethist{ii}]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   modelspec(trainspecs(ii)),[2.4 3.4 3.4 3.0],trainfirsttime,trainrange,.01);

    % now add compaction correction factor
    ms = modelspec(trainspecs(ii));
    ms.thet0 = thetTGG{ii}(1:end-1);
    ms.subfixed = 1:length(ms.thet0);
    [thetTGG{ii},trainsubsubset{ii},logp(ii),thist]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   ms,[3.4 3.0],trainfirsttime(end),trainrange(end),1e6);   
    thethist{ii}=[thethist{ii}; thist];
    
    fid=fopen('thetTGG.tsv','w');
    fprintf(fid,'set\ttraining data\tmodel\tlogp\tN\n');
    for iii=1:length(thetTGG)
        fprintf(fid,[trainlabels{iii} '\t' datasets{trainsets(iii)}.label '\t' ...
                     modelspec(trainspecs(iii)).label]);
        fprintf(fid,['\t(%0.2f)'],logp(iii));
        fprintf(fid,'\t%0.0f',length(trainsubsubset{ii}));
        fprintf(fid,'\t%0.3f',thetTGG{iii});
        fprintf(fid,'\n');
    end
    fclose(fid);
end

% cross-validation
Nk=20;
trainsub=trainsubsubset{1};
trainsubr=trainsub(randperm(length(trainsub)));
ntodo = length(trainsub)/Nk;
partitions=reshape(trainsubr(1:Nk*floor(ntodo)),Nk,[]);
alwaysdo = trainsubr((Nk*floor(ntodo))+1:end);

for ii=1:length(trainspecs)
    wdataset=datasets{trainsets(ii)};
    collinear=modelspec(trainspecs(ii)).subamplinear(1);


    u=unique(wdataset.datid);
    clear fp sdp;
    clear testsitedefp;
    subps=[];
    for pp=1:length(u)
        subp=find(wdataset.datid==u(pp));
        subq=find(wdataset.siteid==u(pp));
        if length(subp)>0
            testtsp{pp}=wdataset.meantime(subp);
            testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                                wdataset.sitecoords(subq,:)];
            testsitedefp.names(pp)=wdataset.sitenames(subq);
            testsitedefp.names2(pp)=wdataset.sitenames(subq);
            testsitedefp.firstage(pp)=min(wdataset.meantime(subp));
            testsitedefp.GISfp(pp)=wdataset.siteGISfp(subq);
            testsitedefp.GIA(pp)=wdataset.siteGIA(subq);
        end
        subps=[subps ; subp];
    end
    [fp(subps),sdp(subps),~,testlocp,logpa,passderivs]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(ii)),thetTGG{ii},trainsub,ones(1,length(thetTGG{ii})),testtsp,refyear,collinear);
    
    clear fp1 sdp1;
    clear testsitedefp1;
    fp1 = zeros(size(fp));
    sdp1 = zeros(size(sdp));
    
    clear subps omit trainsubdo;
    for kk=1:Nk
        trainsubdo{kk} = partitions(setdiff(1:Nk,kk),:);
        trainsubdo{kk} = sort(union(trainsubdo{kk}(:)',alwaysdo));
        omit{kk} = sort(partitions(kk,:));

        u = unique(wdataset.datid(omit{kk}));
        subps{kk}=[];
        
        for pp=1:length(u)
            subp=intersect(find(wdataset.datid==u(pp)),omit{kk});
            subq=find(wdataset.siteid==u(pp));
            if length(subp)>0
                testsitedefp1{kk}.testtsp1{pp}=wdataset.meantime(subp);
                testsitedefp1{kk}.sites(pp,:)=[wdataset.siteid(subq) ...
                                    wdataset.sitecoords(subq,:)];
                testsitedefp1{kk}.names(pp)=wdataset.sitenames(subq);
                testsitedefp1{kk}.names2(pp)=wdataset.sitenames(subq);
                testsitedefp1{kk}.firstage(pp)=min(wdataset.meantime(subp));
                testsitedefp1{kk}.GISfp(pp)=wdataset.siteGISfp(subq);
                testsitedefp1{kk}.GIA(pp)=wdataset.siteGIA(subq);
            end
            subps{kk}=[subps{kk} ; subp];
        end
        
        
    end
    
    clear logpps;
    for kk=1:Nk
        disp(kk);
        
        [intr,ij,ik] = intersect(trainsub,trainsubdo{kk});
        pd=passderivs;
        pd.dK=pd.dK(ij);
        pd.yoffset=pd.yoffset(ij);
        pd.df=pd.df(ij);
        
        [fp1(subps{kk}),sdp1(subps{kk}),Vp1,testlocp1]=RegressHoloceneDataSets(wdataset,testsitedefp1{kk},modelspec(trainspecs(ii)),thetTGG{ii},trainsubdo{kk},ones(1,length(thetTGG{ii})),testsitedefp1{kk}.testtsp1,refyear,collinear,pd);
        
        [alfa,~,s,r]=svdinv(Vp1,wdataset.Y(subps{kk})-fp1(subps{kk})');          
        logpps(kk) = -.5 * (wdataset.Y(subps{kk})'-fp1(subps{kk}))*alfa - .5 * log(2*pi) * length(subps{kk}) - .5 * sum(log(s(1:r)));
        
    end
    logpp(ii)=sum(logpps);
    
    fid=fopen('modelCV.tsv','w');
    fprintf(fid,'set\ttraining data\tmodel\tlogp\tN\tNk\tlog posterior predictive\n');
    for iii=1:length(trainspecs)
        fprintf(fid,[trainlabels{iii} '\t' datasets{trainsets(iii)}.label '\t' ...
                     modelspec(trainspecs(iii)).label]);
        fprintf(fid,['\t(%0.2f)'],logp(iii));
        fprintf(fid,'\t%0.0f',length(trainsubsubset{ii}));
        fprintf(fid,'\t%0.0f',Nk);
        fprintf(fid,'\t%0.0f',logpp(iii));
        fprintf(fid,'\n');
    end
    fclose(fid);
end


save thetTGG thetTGG trainsubsubset


%% now do a regression
clear Loc cnt oldest testsitedef;
preserveall={'North Carolina','New Jersey','Florida','Nova Scotia'};
for ii=1:length(PX.sitenames)
    s0=0;
    for jj=1:length(preserveall)
        s0=s0+length(strfind(PX.sitenames{ii},preserveall{jj}));
    end
    s=strfind(PX.sitenames{ii},'-');
    if (s>0).*(s0==0)
        Loc{ii}=PX.sitenames{ii}(1:s-1);
    else
        Loc{ii}=PX.sitenames{ii};
    end
    sub=strfind(Loc{ii},' ');
    Loc{ii}=Loc{ii}(setdiff(1:length(Loc{ii}),sub));
    sub=find(PX.datid==PX.siteid(ii));
    cnt(ii)=length(sub);
    if length(sub)>0
        oldest(ii)=min(union(PX.time1(sub),PX.time2(sub)));
    else
        oldest(ii)=2000;
    end
end
[uLoc,ui]=unique(Loc);

clear testsitedef;
testsitedef.sites=[0 1e6 1e6];
testsitedef.names={'GSL'};
testsitedef.names2={'GSL'};
testsitedef.firstage=min(oldest);

for ii=1:length(uLoc)
    sub=find(strcmpi(uLoc{ii},Loc));
    [m,mi]=max(cnt(sub));
    if m>0
        si=find(PX.datid==PX.siteid(sub(mi))); si=si(1);
        testsitedef.sites(end+1,:)=[PX.datid(si) PX.lat(si) PX.long(si)];
        testsitedef.names2={testsitedef.names2{:}, PX.sitenames{sub(mi)}};
        testsitedef.names={testsitedef.names{:}, uLoc{ii}};
        testsitedef.firstage = [testsitedef.firstage min(oldest(sub))];
    end
end
for si=1:length(TGNOCW.siteid)
        testsitedef.sites(end+1,:)=[TGNOCW.siteid(si) TGNOCW.sitecoords(si,:)];
        testsitedef.names2={testsitedef.names2{:}, TGNOCW.sitenames{si}};
        u=TGNOCW.sitenames{si};
        testsitedef.names={testsitedef.names{:}, u(setdiff(1:length(u),strfind(u,' ')))};
        
        sub=find(TGNOCW.datid==TGNOCW.siteid(si));
        TGoldest=min(TGNOCW.meantime(sub));
        testsitedef.firstage = [testsitedef.firstage TGoldest];
end


GISfpt.lat=GISfplat;
GISfpt.long=GISfplong;
GISfpt.fp=GISfp;

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;

for ii=1:length(testsitedef.sites(:,1))
    testsitedef.GISfp = interp2(GISfpt.long,GISfpt.lat,GISfpt.fp,testsitedef.sites(:,3),testsitedef.sites(:,2),'linear');
    testsitedef.GISfp(find(testsitedef.sites(:,2)>100))=1;

    testsitedef.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedef.sites(:,2),testsitedef.sites(:,3),'linear');
    testsitedef.GIA(find(testsitedef.sites(:,2)>100))=0;

end


testt = [-1000:20:2000 2010];

regresssets=[1 2 5 1 5 7];
regressparams=[1 1 1 2 2 3];
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);

    noiseMasks = ones(6,length(thetTGG{trainspecs(jj)}));
    noiseMasks(2,[modelspec(trainspecs(jj)).subamplinear modelspec(trainspecs(jj)).subampoffset]  )=0; %without linear
    noiseMasks(3,[modelspec(trainspecs(jj)).subampglobal modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local
    noiseMasks(4,[modelspec(trainspecs(jj)).subamplinear modelspec(trainspecs(jj)).subampglobal  modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local non-linear
    noiseMasks(5,[modelspec(trainspecs(jj)).subampglobal modelspec(trainspecs(jj)).subampnonlinear modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local linear
    noiseMasks(6,setdiff(modelspec(trainspecs(jj)).subamp,modelspec(trainspecs(jj)).subampGIS))=0; %Greenland only
    noiseMasklabels={'full','nonlin','regloc','reglocnonlin','regloclin','gis'};

    wdataset=datasets{ii};

    labls{iii}=['_' regresslabels{iii}];
    disp(labls{iii});

    trainsub = find((wdataset.limiting==0)); % only index points
    wdataset.dY = sqrt(datasets{ii}.dY.^2 + (thetTGG{jj}(end)*wdataset.compactcorr).^2);
    wdataset.Ycv = datasets{ii}.Ycv + diag(thetTGG{jj}(end)*wdataset.compactcorr).^2;
    subtimes=find(testt>=min(union(wdataset.time1,wdataset.time2)));
    
    collinear=modelspec(trainspecs(jj)).subamplinear(1);
    [f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},logp(ii,jj),passderivs,invcv]=RegressHoloceneDataSets(wdataset,testsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);

    labl=labls{iii}; disp(labl);
    
    u=unique(wdataset.datid);
    clear fp sdp;
    clear testsitedefp;
    subps=[];
    for pp=1:length(u)
        subp=find(wdataset.datid==u(pp));
        subq=find(wdataset.siteid==u(pp));
        if length(subp)>0
            testtsp{pp}=wdataset.meantime(subp);
            testsitedefp.sites(pp,:)=[wdataset.siteid(subq) ...
                                wdataset.sitecoords(subq,:)];
            testsitedefp.names(pp)=wdataset.sitenames(subq);
            testsitedefp.names2(pp)=wdataset.sitenames(subq);
            testsitedefp.firstage(pp)=min(wdataset.meantime(subp));
            testsitedefp.GISfp(pp)=wdataset.siteGISfp(subq);
            testsitedefp.GIA(pp)=wdataset.siteGIA(subq);
        end
        subps=[subps ; subp];
    end
    [fp(subps),sdp(subps),~,testlocp]=RegressHoloceneDataSets(wdataset,testsitedefp,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);

    fid=fopen(['TGandProxyData' labl '.tsv'],'w');
    fprintf(fid,['Site\tID\ttime1\ttime2\tlimiting\tGIAproj\tY-GIAproj\tY\' ...
                 'tdY\tcompactcorr\tistg\tlat\tlong\tsite GIA\tmean time\tYposterior\tdYposterior\n']);
    for i=1:size(wdataset.datid,1)
        subq=find(wdataset.siteid==wdataset.datid(i));
        compactcorr=full(wdataset.compactcorr);
        subq=find(wdataset.siteid==wdataset.datid(i));
        fprintf(fid,[wdataset.sitenames{subq} '\t']);
        fprintf(fid,'%d\t',wdataset.datid(i));
        fprintf(fid,'%0.1f\t',wdataset.time1(i));
        fprintf(fid,'%0.1f\t',wdataset.time2(i));
        fprintf(fid,'%d\t',wdataset.limiting(i));
        fprintf(fid,'%0.2f\t',wdataset.GIAproj(i));
        fprintf(fid,'%0.2f\t',wdataset.Y(i));
        fprintf(fid,'%0.2f\t',wdataset.Y0(i));
        fprintf(fid,'%0.2f\t',wdataset.dY(i));
        fprintf(fid,'%0.2f\t',compactcorr(i));
        fprintf(fid,'%d\t',wdataset.istg(i));
        fprintf(fid,'%0.2f\t',wdataset.lat(i));
        fprintf(fid,'%0.2f\t',wdataset.long(i));
        fprintf(fid,'%0.2f\t',wdataset.siteGIA(subq(1)));
        fprintf(fid,'%0.1f\t',wdataset.meantime(i));
        fprintf(fid,'%0.2f\t',fp(i));
        fprintf(fid,'%0.2f\n',sdp(i));
    end
    fclose(fid)    
    
    if dosldecomp
        makeplots_sldecomp(wdataset,f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},labl);
    end
    
    testreg=testlocs{ii,jj}.reg;
    testsites=testlocs{ii,jj}.sites;
    testX=testlocs{ii,jj}.X;
    testnames2=testlocs{ii,jj}.names2;

    firstyears=[ 0 1000  200  -500  -300 700 1700 1800 1900 900 900  1000 1000 1500 1500 860  1000 1080 1460 1560 1800 1860];
    lastyears=[1800 1800 1000 1000 700 1700 1800 1900 2000 1400 1500 1400 1500 1800 1860 1560 1400 1360 1880 1800 1880 1900];

     [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites,testreg,testX(:,3),firstyears,lastyears);
     [fslopelin,sdslopelin]=SLRateCompare(f2s{ii,jj}(:,5),V2s{ii,jj}(:,:,5),testsites,testreg,testX(:,3),testt(end-1),testt(end));
    
    fid=fopen(['linrates' labl '.tsv'],'w');
    fprintf(fid,['Regional+Local Linear Rates (mm/y), ' labl '\n']);
    fprintf(fid,'Site\tSiteID\tLat\tLong\tICE5G VM2-90\tRate (linear)\t2s');
    for pp=1:length(firstyears)
        fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
    end
    for pp=1:length(diffplus)
        fprintf(fid,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
                     '%0.0f-%0.0f)\t2s'],[firstyears(diffplus(pp)) ...
                            lastyears(diffplus(pp)) firstyears(diffless(pp)) lastyears(diffless(pp))]);
    end
    
    fprintf(fid,'\n');
    for kk=1:size(testsites,1)
        fprintf(fid,testnames2{kk});
        fprintf(fid,'\t%0.2f',testsitedef.sites(kk,:));
        fprintf(fid,'\t%0.2f',testsitedef.GIA(kk));
        fprintf(fid,'\t%0.2f',[fslopelin(kk) 2*sdslopelin(kk)]);
        for pp=1:length(firstyears)
            fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) 2*sdslopeavg(kk,pp)]);
        end
        for pp=1:length(diffplus)
            fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
        end
        fprintf(fid,'\n');
    end
    fclose(fid);

    %%%%%
    
    firstyears=[ 0  1000  200 -300 700 1700 1800 1900 ];
    lastyears=[1800 1800 1000 700 1700 1800 1900 2000 ];

    for pp=1:length(firstyears)
        
        [fslope]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites,testreg,testX(:,3),firstyears(pp),lastyears(pp));
        labl2=[num2str(firstyears(pp)) '-' num2str(lastyears(pp))];
        
        figure;
        worldmap('world');
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        geoshow('landareas.shp','facecolor','none');
  
        %
        %        plotcont; hold on;

         % plot tide gauge locations
        sub=find((testsites(:,2)<1e3).*(testsites(:,1)<5e3));
        sub1=intersect(sub,find(isnan(fslope)));
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),10,[.5 .5 .5],'filled','d'); hold on;
        %        hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),10,[.5 .5 .5],'d'); hold on;
        hq=getkids(hq);
        set(hq,'linew',1);
        sub1=intersect(sub,find(~isnan(fslope)));
        %        scatter(mod(testsites(sub1,3),360),testsites(sub1,2),10,fslope(sub1),'filled','d'); hold on;
        scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),10,fslope(sub1),'filled','d'); hold on;
        
        % plot proxy sites
        sub=find((testsites(:,2)<1e3).*(testsites(:,1)>=5e3).*(testsites(:,1)<1e6));
        sub1=intersect(sub,find(isnan(fslope)));
        %hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),20,[.5 .5 .5]); hold on;
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),20,[.5 .5 .5],'filled'); hold on;
        hq=getkids(hq);
        set(hq,'linew',1);
        sub1=intersect(sub,find(~isnan(fslope)));
        %scatter(mod(testsites(sub1,3),360),testsites(sub1,2),30,fslope(sub1),'filled'); hold on;
        scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),30,fslope(sub1),'filled'); hold on;

        % plot proxy sites (EH12)
        sub=find((testsites(:,1)>=1e6));        
        sub1=intersect(sub,find(~isnan(fslope)));
        %hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),10,fslope(sub1),'filled','s'); hold on;
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),10,fslope(sub1),'filled','s'); hold on;
        
        axis tight;
        colorbar;
        
        %        box on;
        if firstyears(pp)<1900
           caxis([-0.5 2.0]);
        end
        
        title(['Sea level  rates, ' labl2 ' CE (mm/y)']);
        pdfwrite(['map_linrates_global' labl '_' labl2]);

        
        %%%%
        clf;
        worldmap([27 49],[265 300]);
        setm(gca,'parallellabel','off','meridianlabel','off','flinewidth',1);
        ax=gca;
        %geoshow('landareas.shp','facecolor','none');
        usstates = shaperead('usastatehi', 'UseGeoCoords', true);
        stateColor='none';
        geoshow(ax(1), usstates,  'FaceColor', stateColor);
        if exist([IFILES '/province/PROVINCE.SHP']); canada = shaperead([IFILES '/province/PROVINCE.SHP'],'UseGeoCoords',true); else; canada=[]; end
        geoshow(ax(1), canada,  'FaceColor', stateColor);
        setm(ax, 'ParallelLabel', 'on', 'MeridianLabel', 'on')


        %figure;
        %usstates = shaperead('usastatehi', 'UseGeoCoords', true);
        %plot((mod([usstates.Lon],360)),[usstates.Lat],'k','linew',.5); hold on;
        %if exist([IFILES '/province/PROVINCE.SHP']); canada = shaperead([IFILES '/province/PROVINCE.SHP'],'UseGeoCoords',true); else; canada=[]; end
        %if length(canada)>0 ; plot((mod([canada.Lon],360)),[canada.Lat],'k','linew',.5); end
        %xl=xlim;
        %yl=ylim;

        % plot tide gauge locations
        sub=find((testsites(:,2)<1e3).*(testsites(:,1)<5e3));
        sub1=intersect(sub,find(isnan(fslope)));
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),30,[.5 .5 .5],'filled','d'); hold on;
        %        hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),30,[.5 .5 .5],'d'); hold on;
        hq=getkids(hq);
        set(hq,'linew',1);
        sub1=intersect(sub,find(~isnan(fslope)));
        %        scatter(mod(testsites(sub1,3),360),testsites(sub1,2),30,fslope(sub1),'filled','d'); hold on;
        scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),30,fslope(sub1),'filled','d'); hold on;
        
        % plot proxy sites
        sub=find((testsites(:,2)<1e3).*(testsites(:,1)>=5e3).*(testsites(:,1)<1e6));
        sub1=intersect(sub,find(isnan(fslope)));
        %hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),70,[.5 .5 .5]); hold on;
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),70,[.5 .5 .5],'filled'); hold on;
        hq=getkids(hq);
        set(hq,'linew',1);
        sub1=intersect(sub,find(~isnan(fslope)));
        %scatter(mod(testsites(sub1,3),360),testsites(sub1,2),100,fslope(sub1),'filled'); hold on;
        scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),100,fslope(sub1),'filled'); hold on;

        % plot proxy sites (EH12)
        sub=find((testsites(:,1)>=1e6));        
        sub1=intersect(sub,find(~isnan(fslope)));
        %hq=scatter(mod(testsites(sub1,3),360),testsites(sub1,2),30,fslope(sub1),'filled','s'); hold on;
        hq=scatterm(testsites(sub1,2),mod(testsites(sub1,3),360),30,fslope(sub1),'filled','s'); hold on;
        
        axis tight;
        colorbar;
        if firstyears(pp)<1900
            caxis([-0.5 2]);
        end

        %box on;
        title(['Sea level  rates, ' labl2 ' CE (mm/y)']);
        
        %       xlim([265 300]);
        %ylim([27 49]);

        pdfwrite(['map_linrates_NA' labl '_' labl2]);
        
    end

    %%% 

    refyear=2000;
    Mref = eye(size(testX,1));
    for i=1:size(testsites,1)
        
        sub1=find(testreg==testsites(i,1));
        sub2=intersect(sub1,find(testX(:,3)==refyear));
        
        Mref(sub1,sub2)=Mref(sub1,sub2)-1;

    end
    Mref=sparse(Mref);
    
    refyear2=1900;
    Mref2 = eye(size(testX,1));
    for i=1:size(testsites,1)
        
        sub1=find(testreg==testsites(i,1));
        sub2=intersect(sub1,find(testX(:,3)==refyear2));
        
        Mref2(sub1,sub2)=Mref2(sub1,sub2)-1;

    end
    Mref2=sparse(Mref2);    
    
    selsitenames={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','NovaScotia-Chezzetcook'};
    shortnames = {'FL','NC','NJ','NS'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end
    colrs={'r','b','m','g'};

    datsub=find(ismember(testreg,testsites(sitesub,1)));

    % figure of sites with data, full curves
    % figure of sites with data, detrended

    for selmask=[1]
        for dodetrend=[0 1]
            clf;
            
            wf=Mref*f2s{ii,jj}(:,selmask);
            wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
            wsd=sqrt(diag(wV));
            
            labl2='';
            
            if dodetrend
                [wf,wV,wsd]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[0 1000 1700],1800,refyear);
                labl2='_detrended';
            end
            
            [hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,200,testnames2(sitesub));
            
            if selmask==1
                legend(hl,shortnames,'Location','Southwest');
            else
                legend(hl,shortnames,'Location','Southwest');
            end
            set(hp,'xlim',[-500 2010]);
            set(hp(2),'ylim',[-.5 3.2]);
            delete(hp(2));
            

            pdfwrite(['paleorsl_' noiseMasklabels{selmask} labl labl2]);

            fid=fopen(['paleorsl_' noiseMasklabels{selmask} labl labl2 '.tsv'],'w');
            fprintf(fid,outtable);
            fclose(fid);
        end
        
    end

    % figure of GSL and rate of GSL change

    selsitenames={'GSL'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end
    colrs={'k'};

    selmask=1;

    wf=Mref*f2s{ii,jj}(:,selmask);
    wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
    wsd=sqrt(diag(wV));

    wf1900=Mref2*f2s{ii,jj}(:,selmask);
    wV1900=Mref2*V2s{ii,jj}(:,:,selmask)*Mref2';
    wsd1900=sqrt(diag(wV1900));


    
    for dodetrend=[0 1]
        labl3='';
        if dodetrend
            [wf,wV,wsd]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[1000],1800,refyear);
            labl3='_detrended';
        end
        
    
        datsub=find(ismember(testreg,testsites(sitesub,1)));

        for timesteps=[100 400 60 40 20]

            clf;
            [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,timesteps,{'GSL'});
            set(hp,'xlim',[-500 2010]);
            
            labl2=[labl labl3];  

            delete(hp(2));
            if timesteps==100
                pdfwrite(['GSL_' num2str(timesteps) 'y' labl2]);
            end
            


            fid=fopen(['GSL_' num2str(timesteps) 'y' labl2 '.tsv'],'w');
            fprintf(fid,outtable);

            Nsamps=10000;
            sub1=find((difftimes<=1800).*(difftimes>=1500)); sub2=find(difftimes>1800);
            sub1a=find((difftimes<=1800).*(difftimes>=0)); sub2=find(difftimes>1800); 
            if ((length(sub1a))>0).*(length(sub1)>0)
                samps=mvnrnd(dGSL,dGSLV,Nsamps);
                if (length(sub1)>0).*(length(sub2)>0)
                    [max1,max1i]=max(samps(:,sub1),[],2);
                    [max1a,max1ia]=max(samps(:,sub1a),[],2);
                    difr=bsxfun(@minus,samps(:,sub2),max1);
                    difra=bsxfun(@minus,samps(:,sub2),max1a);
                    Ppositive = sum(samps(:,sub2)>0,1)/size(samps,1);
                    q=cumprod((samps(:,sub2(end:-1:1))>0)*1,2); q=q(:,end:-1:1);
                    Ppositiveseries = sum(q,1)/size(samps,1);
                    Plarger = sum(difr>0,1)/size(samps,1);
                    Plargera = sum(difra>0,1)/size(samps,1);
                    fprintf(fid,'\n\nProbability faster than during fastest interval\n');
                    fprintf(fid,'Central year\tP > 0\tP all subsequent > 0\tP centered 1500-1800\tP centered 0-1800');
                    fprintf(fid,'\n%0.0f\t%0.3f\t%0.3f\t%0.3f\t%0.3f',[difftimes(sub2)' ; Ppositive; Ppositiveseries ; Plarger ; Plargera]);
                end
            end
            
            suboverallyrs=find(mod(difftimes,timesteps)==timesteps/2);
            subyrs=find(difftimes(suboverallyrs)>=1800);
            wquants=[.67 .9 .95 .99 .995 .999];
            if length(subyrs)>0
                fprintf(fid,'\n\nLast year faster');
                fprintf(fid,'Central year');
                fprintf(fid,'\t%0.2f',wquants);
                fprintf(fid,'\n');
                
                clear lastyrlarger;
                for mm=1:length(subyrs)
                    curyr=difftimes(suboverallyrs(subyrs(mm)));
                    for nn=1:Nsamps
                        sub=intersect(suboverallyrs,find(difftimes<curyr));
                        sub2=find(samps(nn,sub)>samps(nn,suboverallyrs(subyrs(mm))));
                        if length(sub2)>0
                            lastyrlarger(nn,mm)=difftimes(suboverallyrs(sub2(end)));
                        else
                            lastyrlarger(nn,mm)=min(difftimes)-1;
                        end
                    end
                    fprintf(fid,'\n%0.0f',curyr);
                    fprintf(fid,'\t%0.0f',quantile(lastyrlarger(:,mm),wquants));
                end
            end
            fclose(fid);
            
     
        end
           

        % output GSL covariance
        

        fid=fopen(['GSL'  labl2 '_cov.tsv'],'w');
        fprintf(fid,'mm^2');
        fprintf(fid,'\t%0.0f',testX(datsub,3));
        fprintf(fid,'\n');
        for ppp=1:length(datsub)
            fprintf(fid,'%0.0f',testX(datsub(ppp),3));
            fprintf(fid,'\t%0.3e',wV(datsub(ppp),datsub));
            fprintf(fid,'\n');
        end
        
        fclose(fid);
    end

    
    [~,~,~,~,~,~,outtable1900]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf1900(datsub),wV1900(datsub,datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,100,{'GSL'});
    
    labl2=labl;
    fid=fopen(['GSL_' num2str(timesteps) 'y' labl2 '_ref1900.tsv'],'w');
    fprintf(fid,outtable1900);
    fclose(fid);
   
    
    %%%

    % figure of GIS and rate of GIS change

    selsitenames={'GSL'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end
    colrs={'k','r','c','b','g','m','y'};

    selmask=6;

    wf=Mref*f2s{ii,jj}(:,selmask);
    wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
    wsd=sqrt(diag(wV));

    datsub=find(ismember(testreg,testsites(sitesub,1)));

    for timesteps=[100 600]

        clf;
        [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,timesteps,{'GIS'});
        set(hp,'xlim',[-1000 2010]);

        pdfwrite(['GIS_' num2str(timesteps) 'y' labl]);


        fid=fopen(['GIS_' num2str(timesteps) 'y' labl '.tsv'],'w');
        fprintf(fid,outtable);
        
    end

    %%%%

    % figure of sea-level differences

    selsitenames={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','NovaScotia-Chezzetcook','NEWYORK','HALIFAX','Connecticut'};
    shortselnames={'FL','NC','NJ','NS','NYC','HFX','CT'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end
    colrs={'r','m','g','b','c','k','y'};

    selmask=1;

    clear pairsets;
    pairsets{1}=[1 2; 2 3; 2 4; 3 4; 5 3; 5 4 ; 7 3];
    pairsets{2}=[1 2; 1 3 ; 1 4 ; 2 3; 2 4 ; 3 4];
    pairsets{3}=[1 2 ; 3 2 ; 4 3]
    for mmm=1:length(pairsets)
        for dodetrend=[0 1]
            
            wdiffpair=[]; wpairnames={};
            gradf=[]; gradsd=[]; gradt=[]; gradpair=[]; gradV=[]; ...
                  gradstarttimes=[]; clear pairnames;
            diffpairs=pairsets{mmm};
            labl2=['_' num2str(mmm)];
            wf=f2s{ii,jj}(:,selmask);
            wV=V2s{ii,jj}(:,:,selmask);

            if dodetrend
                [wf1,wV1,wsd1,sitespec]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[0],1800,refyear);
                [wf2,wV2,wsd2]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[1000],1800,refyear);
                
                labl2=[labl2 '_detrended'];
            end

            for i=1:size(diffpairs,1)
                pairnames{i} = [shortselnames{diffpairs(i,2)} '-' shortselnames{diffpairs(i,1)}];
                sub1=find(testreg==testsites(sitesub(diffpairs(i,1))));
                sub2=find(testreg==testsites(sitesub(diffpairs(i,2))));
                [u,ui,uj]=intersect(testX(sub1,3),testX(sub2,3));
                Mgrad=sparse(length(u),length(testX));
                Mgrad(:,sub1(ui))=-eye(length(ui));
                Mgrad(:,sub1(ui(end))) = Mgrad(:,sub1(ui(end)))+1;
                
                Mgrad(:,sub2(uj)) = Mgrad(:,sub2(uj)) + eye(length(uj));
                Mgrad(:,sub2(uj(end))) = Mgrad(:,sub2(uj(end)))-1;
                
                if dodetrend
                    %                   wf1=wf2; wV1=wV2; labl2=['_' num2str(mmm) '_detrended1000']; ; % always detrend with respect to 1000-1800
                    wf = wf1; wV = wV1;
                    q = Mgrad*wf;
                    if sum(isnan(q))
                        wf=wf2; wV=wV2;
                        q=Mgrad*wf;
                    end
                else
                    q=Mgrad*wf;
                end
                
                
                if sum(~isnan(q))>1

                    gradf = [gradf ; q];
                    gradV(length(gradV) + [1:length(u)],length(gradV) + [1:length(u)]) = Mgrad*wV*Mgrad';
                    gradt = [gradt ; u];
                    gradpair = [gradpair ; ones(length(u),1)*i];
                    gradstarttimes=[gradstarttimes ; u(1)];
                    wdiffpair = [wdiffpair ; i];
                    wpairnames={wpairnames{:},pairnames{i}};
                end
                
                
                
            end 
            
            gradsd = sqrt(diag(gradV));

            %       [gradf,gradV,gradsd]=DetrendSLReconstruction(gradf,gradV,[1:length(diffpairs)]',gradpair,gradt,[0 1000],1800,refyear);
            
            for timesteps=[100 60 40 20 ]

                clf;
                [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(gradt,gradpair,wdiffpair,gradf,gradV,colrs,gradstarttimes,2010,0,timesteps,wpairnames);
                set(hp,'xlim',[-500 2000]);
                if mmm>1
                    set(hp,'xlim',[0 2010]);
                end
                legend(hl,wpairnames,'Location','Southwest');
                if dodetrend
                    delete(hp(2));
                end
                
                if timesteps==100
                    pdfwrite(['RSLgrad_' num2str(timesteps) 'y' labl labl2]);
                end
                
                
                
                %if (mmm==1)
                    firstyears=[ 0   1000  200 700  1700 1800 1900 900  1000 1600 1860 1940 1980];
                    lastyears= [1800 1800 1000 1700 1800 1900 2000 1500 1300 1800 2000 1980 2010];
                    [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(gradf,gradV,wdiffpair,gradpair,gradt,firstyears,lastyears);
                    
                    fid=fopen(['RSLgrad_' num2str(timesteps) 'y' labl labl2 '.tsv'],'w');
                    
                    
                    fprintf(fid,'Pair');
                    for pp=1:length(firstyears)
                        fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
                    end
                    for pp=1:length(diffplus)
                        fprintf(fid,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
                                     '%0.0f-%0.0f)\t2s'],[firstyears(diffplus(pp)) ...
                                            lastyears(diffplus(pp)) firstyears(diffless(pp)) lastyears(diffless(pp))]);
                    end
                    
                    fprintf(fid,'\n');
                    for kk=1:length(wpairnames)
                        fprintf(fid,wpairnames{kk});
                        for pp=1:length(firstyears)
                            fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) 2*sdslopeavg(kk,pp)]);
                        end
                        for pp=1:length(diffplus)
                            fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
                        end
                        fprintf(fid,'\n');
                    end
                    fprintf(fid,'\n\n');
                    
                    fprintf(fid,outtable);
                    fclose(fid);
                    % end
                
                
                
            end
            
        end
    end
    
        % weight GIA models
        
        GIAtimes=[.15 .95 1.95];
        GIAtimes=union(GIAtimes,[0 .05:.2:3.05 ]);
        [GIAt,GIAsl,GIAsites,GIAicehist,GIAsolidearth] = readJXMSiteGIA(GIAfiles,GIAtimes);
        
        t0 = find(GIAt==.15);
        t1 = find(GIAt==1.05);
        t2 = find(GIAt==1.95);
        t3 = find(GIAt==0.25);
        t4 = find(GIAt==0.05);
        
        GIAavg = squeeze([GIAsl(t0,:,:)-GIAsl(t2,:,:)])/-(.15-1.95);       
        GIAavgA = squeeze([GIAsl(t0,:,:)-GIAsl(t1,:,:)])/-(.15-1.05);
        GIAavgB = squeeze([GIAsl(t1,:,:)-GIAsl(t2,:,:)])/-(1.05-1.95);
        GIAavgZ = squeeze([GIAsl(t4,:,:)-GIAsl(t3,:,:)])/(.25-.05);
        
        GIAsitetarg={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','ChristmasIsland-Multiple','Massachusetts-WoodIsland','EH12_1','EH12_2','EH12_3','EH12_4','EH12_5','EH12_6','EH12_7','EH12_8','EH12_9','EH12_10','EH12_11','EH12_13','EH12_14','EH12_16'};
        %GIAsitetarg={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','ChristmasIsland-Multiple','EH12_6','EH12_7','EH12_8','EH12_9','EH12_10','EH12_11','EH12_13','EH12_14','EH12_16'};
        
        clear GIAsub;
        for pp=1:length(GIAsitetarg)
            GIAsub(pp)=find(strcmpi(GIAsitetarg{pp},GIAsites));
        end
        GIAavg=GIAavg(GIAsub,:);
        GIAavgA=GIAavgA(GIAsub,:);
        GIAavgB=GIAavgB(GIAsub,:);
        GIAavgZ=GIAavgZ(GIAsub,:);
        
                selsitenames={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','ChristmasIsland','Massachusetts','EH12_1','EH12_2','EH12_3','EH12_4','EH12_5','EH12_6','EH12_7','EH12_8','EH12_9','EH12_10','EH12_11','EH12_13','EH12_14','EH12_16'};
                %selsitenames={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','ChristmasIsland','EH12_6','EH12_7','EH12_8','EH12_9','EH12_10','EH12_11','EH12_13','EH12_14','EH12_16'};
        sitesub=[];
        for kk=1:length(selsitenames)
            q=find(strcmpi(selsitenames{kk},testsitedef.names));
            if length(q)>0
                sitesub=[sitesub q(1)];
            end
        end
        [fsl,Vsl]=SLRateMultisite(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites(sitesub,:),testreg,testX(:,3),1000,1800);
        [fslA,sdslA,fslAd,sdslAd]=SLRateCompare(f2s{ii,jj}(:,1),V2s{ii,jj}(:,:,1),testsites(sitesub,:),testreg,testX(:,3),[0 900],[900 1800])

        Mdiff=zeros(length(sitesub)-1,length(sitesub)); Mdiff(:,1)=-1; Mdiff(:,2:end)=eye(length(sitesub)-1);
        
        GIAavg2=Mdiff*GIAavg;
        fsl2=Mdiff*fsl;
        Vsl2=Mdiff*Vsl*Mdiff';
        Vsl2inv=pinv(Vsl2);
        dfsl2=bsxfun(@minus,GIAavg2,fsl2);
        
        logp = -diag(dfsl2'*(Vsl2inv/1)*dfsl2);
        logp=logp-max(logp);
        probwts=exp(logp);
        probwts=probwts/sum(probwts);
        
        GIAwtmean = sum(bsxfun(@times,reshape(probwts,1,1,[]),GIAsl),3);
        GIAwtrate = (GIAwtmean(t0,:)-GIAwtmean(t2,:))/-(.15-1.95);
        GIAwtmeandetr = GIAwtmean+bsxfun(@times,GIAwtrate,GIAt'-.15); GIAwtmeandetr=bsxfun(@minus,GIAwtmeandetr,GIAwtmeandetr(t0,:));
       GIAwtrate2 = (GIAwtmeandetr(t0,:)-GIAwtmeandetr(t2,:))/-(.15-1.95);
        
        GIAsitetarg={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','NovaScotia-Chezzetcook'};
        shortnames={'FL','NC','NJ','NS'};
        diffpairs={[2 1],[3 2],[4 2],[4 3]}; 
        colrs={'r','b','m','g'};
        clear GIAsub;
        clf;
        subplot(2,1,1);
        for pp=1:length(GIAsitetarg)
            GIAsub(pp)=find(strcmpi(GIAsitetarg{pp},GIAsites));
            plot(1950-1000*GIAt,1000*GIAwtmeandetr(:,GIAsub(pp)),colrs{pp},'linew',2); hold on;
        end
        legend(shortnames,'location','southeast');
        xlim([-500 2000]);
        xlabel('Year (CE)'); ylabel('mm'); title('Detrended mean GIA estimate');
        pdfwrite(['GIAwtmeandetr' labl]);

        clf;
        subplot(2,1,1);
        for pp=1:length(GIAsitetarg)
            GIAsub(pp)=find(strcmpi(GIAsitetarg{pp},GIAsites));
            plot(1950-1000*GIAt,1000*GIAwtmean(:,GIAsub(pp)),colrs{pp},'linew',2); hold on;
        end
        legend(shortnames,'location','southeast');
        xlabel('Year (CE)'); ylabel('mm'); title('Mean GIA estimate');
        xlim([-500 2000]);
        pdfwrite(['GIAwtmean' labl]);

        
        
        %%
        
        for dodetrend=[0 1]
            wf=f2s{ii,jj}(:,1); wV=V2s{ii,jj}(:,:,1);
            colrs={'r','b','m','g'};

            sitesub2=[];
            Mdiff=speye(length(wf));
            for kk=1:length(GIAsitetarg)
                q=find(strcmpi(GIAsitetarg{kk},testsitedef.names));
                if length(q)>0
                    sitesub2=[sitesub2 q(1)];
                    datsub=find(testreg==testsitedef.sites(q(1)));
                    wf(datsub,:)=wf(datsub,:)-interp1(1950-GIAt*1000,1000*GIAwtmean(:,GIAsub(kk)),testX(datsub,3),'linear','extrap');
                    basesub=intersect(datsub,find((testX(:,3)>=2000).*(testX(:,3)<=2000)));
                    Mdiff(datsub,basesub)= Mdiff(datsub,basesub)-1/length(basesub);
                end
            end
            
            
            datsub=find(ismember(testreg,testsites(sitesub2)));
            datsub=intersect(datsub,find(~isnan(wf)));

            labl2='';
            
            if dodetrend
                [wf(datsub,:),wV(datsub,datsub)]=DetrendSLReconstruction(wf(datsub,:),wV(datsub,datsub),testsites(sitesub2,:),testreg(datsub),testX(datsub,3),1000,1800,refyear);
                labl2='_detrended';
            end

            wf=Mdiff*wf;
            wV=Mdiff*wV*Mdiff';    

            gradf=[]; gradsd=[]; gradt=[]; gradpair=[]; gradV=[]; ...
                  gradstarttimes=[]; wdiffpair=[]; wpairnames={}; clear pairnames;
            for i=1:length(diffpairs)
                pairnames{i} = [shortnames{diffpairs{i}(2)} '-' shortnames{diffpairs{i}(1)}];
                q1=find(strcmpi(GIAsitetarg{diffpairs{i}(1)},testsitedef.names));
                q2=find(strcmpi(GIAsitetarg{diffpairs{i}(2)},testsitedef.names));
                
                sub1=find(testreg==testsitedef.sites(q1(1)));
                sub2=find(testreg==testsitedef.sites(q2(1)));
                
                [u,ui,uj]=intersect(testX(sub1,3),testX(sub2,3));
                Mgrad=sparse(length(u),length(testX));
                Mgrad(:,sub1(ui))=-eye(length(ui));
                Mgrad(:,sub1(ui(end))) = Mgrad(:,sub1(ui(end)))+1;
                
                Mgrad(:,sub2(uj)) = Mgrad(:,sub2(uj)) + eye(length(uj));
                Mgrad(:,sub2(uj(end))) = Mgrad(:,sub2(uj(end)))-1;
                
                q = Mgrad*wf;
                
                if sum(~isnan(q))>1

                    gradf = [gradf ; q];
                    gradV(length(gradV) + [1:length(u)],length(gradV) + [1:length(u)]) = Mgrad*wV*Mgrad';
                    gradt = [gradt ; u];
                    gradpair = [gradpair ; ones(length(u),1)*i];
                    gradstarttimes=[gradstarttimes ; u(1)];
                    wdiffpair = [wdiffpair ; i];
                    wpairnames={wpairnames{:},pairnames{i}};
                end
                
            end 
            
            clf;
            
            [hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub2,1),wf(datsub),wV(datsub,datsub),colrs,testsitedef.firstage(sitesub2),testt(end),0,200,testnames2(sitesub2));
            
            if selmask==1
                legend(hl,shortnames,'Location','Southwest');
            else
                legend(hl,shortnames,'Location','Southwest');
            end
            set(hp,'xlim',[-500 2010]);
            set(hp(2),'ylim',[-.5 3.2]);
            delete(hp(2));
            if dodetrend
                title(hp(1),'RSL less wt mean GIA - detrended');
            else
                title(hp(1),'RSL less wt mean GIA');
            end
            
            pdfwrite(['rsllessGIAwtmean' labl labl2]);
            
            fid=fopen(['rsllessGIAwtmean_' noiseMasklabels{selmask} labl labl2 '.tsv'],'w');
            fprintf(fid,outtable);
            fclose(fid);
            
            gradsd = sqrt(diag(gradV));

            %       [gradf,gradV,gradsd]=DetrendSLReconstruction(gradf,gradV,[1:length(diffpairs)]',gradpair,gradt,[0 1000],1800,refyear);
            
            for timesteps=[100]

                clf;
                [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(gradt,gradpair,wdiffpair,gradf,gradV,colrs,gradstarttimes,2010,0,timesteps,wpairnames);
                set(hp,'xlim',[-500 2000]);
                legend(hl,wpairnames,'Location','Southwest');
                delete(hp(2));                
                pdfwrite(['rsllessGIAwtmeangrad_' num2str(timesteps) 'y' labl labl2]);
                
                fid=fopen(['rsllessGIAwtmeangrad_' noiseMasklabels{selmask} labl labl2 '.tsv'],'w');
                fprintf(fid,outtable);
                fclose(fid);               
                
            end
            
        end
        
        %%
        
        
        clf;
        subplot(2,1,1);
        clear pairnames;
        for pp=1:length(diffpairs)
            plot(1950-1000*GIAt,1000*(GIAwtmeandetr(:,GIAsub(diffpairs{pp}(1)))-GIAwtmeandetr(:,GIAsub(diffpairs{pp}(2)))),colrs{pp},'linew',2); hold on;
            pairnames{pp}=[shortnames{diffpairs{pp}(1)} '-' shortnames{diffpairs{pp}(2)}];
        end
       
        legend(pairnames,'location','southeast');
        xlabel('Year (CE)'); ylabel('mm'); title('Detrended mean GIA estimate');
        pdfwrite(['GIAwtmeandetr_grad' labl]);
        %%%%%
        
        [lsort,lsorti]=sort(testsitedef.sites(sitesub,2));
% $$$ 
% $$$         clf;
% $$$         subplot(1,2,1);
% $$$         clear hp;
% $$$         hp(1)=plot(GIAwtrate(GIAsub(lsorti)),lsort,'k','linew',2); hold on;
% $$$         for nnnn=1:length(selsitenames)
% $$$             hp(2)=plot(fsl(nnnn),testsitedef.sites(sitesub(nnnn),2),'rd'); hold on;
% $$$             plot(fsl(nnnn)+[-1 1]*2*sqrt(Vsl(nnnn,nnnn)),testsitedef.sites(sitesub(nnnn),2)*[1 1],'r-'); hold on;
% $$$         end
% $$$         xlabel('Avg. rate (0-1800 CE), mm/y');
% $$$         legend(hp,'GIA Model','Proxy','Location','Southwest');
% $$$         ylabel('Latitude');
% $$$         ylim([31 45]); xlim([0 2]);
% $$$         
% $$$         pdfwrite('GIAwtratebylat');
        
        %%%%
        [slogp,slogpi]=sort(logp,'descend');
        fid=fopen(['GIAcompare' labl labl2 '.tsv'],'w');
        fprintf(fid,'\tlat\tlong');
        fprintf(fid,'\trate\t2s');
        for nnnn=slogpi(:)'
            fprintf(fid,['\t' GIAicehist{nnnn} '_' GIAsolidearth{nnnn}]);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Relative log L\t\t\t\t');
        fprintf(fid,'\t%0.3e',slogp);
        fprintf(fid,'\n');
        
        fprintf(fid,'Rates\n');
        for nnnn=1:length(selsitenames)
            fprintf(fid,selsitenames{nnnn});
            fprintf(fid,'\t%0.2f',testsitedef.sites(sitesub(nnnn),2:3));
            fprintf(fid,'\t%0.2f',[fsl(nnnn) 2*sqrt(Vsl(nnnn,nnnn))]);
            fprintf(fid,'\t%0.2f',GIAavg(nnnn,slogpi));
            fprintf(fid,'\n');
        end
        
        
        fprintf(fid,['Rate differences rel. to ' selsitenames{1} '\n']);      
        for nnnn=2:length(selsitenames)
            fprintf(fid,selsitenames{nnnn});
            fprintf(fid,'\t%0.2f',testsitedef.sites(sitesub(nnnn),2:3));
            fprintf(fid,'\t%0.2f',[fsl2(nnnn-1) 2*sqrt(Vsl2(nnnn-1,nnnn-1))]);
            fprintf(fid,'\t%0.2f',GIAavg2(nnnn-1,slogpi));
            fprintf(fid,'\n');
        end
        
        fprintf(fid,'Rates difference, 900-1800 vs 0-900 CE \n');
        for nnnn=1:length(selsitenames)
            fprintf(fid,selsitenames{nnnn});
            fprintf(fid,'\t%0.2f',testsitedef.sites(sitesub(nnnn),2:3));
            fprintf(fid,'\t%0.2f',[fslAd(nnnn) 2*sdslAd(nnnn)]);            
            fprintf(fid,'\t%0.2f',GIAavgA(nnnn,slogpi)-GIAavgB(nnnn,slogpi));
            fprintf(fid,'\n');
        end        
        
        fprintf(fid,'Rates, 1700-1900 CE\t\t\t\tICE5G VM2-90 \n');
        for nnnn=1:length(selsitenames)
            fprintf(fid,selsitenames{nnnn});
            fprintf(fid,'\t%0.2f',testsitedef.sites(sitesub(nnnn),2:3));
            fprintf(fid,'\t');            
            fprintf(fid,'\t%0.2f',testsitedef.GIA(sitesub(nnnn)));            
            fprintf(fid,'\t%0.2f',GIAavgZ(nnnn,slogpi));
            fprintf(fid,'\n');
        end        
                
        fclose(fid);
        
        
                
        %%%%
        % Do GIA comparison 
        t0 = find(GIAt==.15);
        t1 = find(GIAt==1.05);
        t2 = find(GIAt==1.95);
        t3 = find(GIAt==0.25);
        t4 = find(GIAt==0.05);
        
        GIAavg2 = squeeze([GIAsl(t0,:,:)-GIAsl(t2,:,:)])/-(.15-1.95);       
        GIAavgA2 = squeeze([GIAsl(t0,:,:)-GIAsl(t1,:,:)])/-(.15-1.05);
        GIAavgB2 = squeeze([GIAsl(t1,:,:)-GIAsl(t2,:,:)])/-(1.05-1.95);
        GIAavgZ2 = squeeze([GIAsl(t4,:,:)-GIAsl(t3,:,:)])/(.25-.05);
        
        
        [slogp,slogpi]=sort(logp,'descend');
        fid=fopen(['GIAmodels' labl labl2 '.tsv'],'w');
        for nnnn=slogpi(:)'
            fprintf(fid,['\t' GIAicehist{nnnn} '_' GIAsolidearth{nnnn}]);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Relative log L');
        fprintf(fid,'\t%0.3e',slogp);
        fprintf(fid,'\n');
        
        fprintf(fid,'Rates 0-1800\n');
        for nnnn=1:length(GIAsites)
            fprintf(fid,GIAsites{nnnn});
            fprintf(fid,'\t%0.2f',GIAavg2(nnnn,slogpi));
            fprintf(fid,'\n');
        end       
        
        fprintf(fid,'\n\nRates 1700-1900 \n');
         for nnnn=1:length(GIAsites)
            fprintf(fid,GIAsites{nnnn});
            fprintf(fid,'\t%0.2f',GIAavgZ2(nnnn,slogpi));
            fprintf(fid,'\n');
         end
         
                
        fclose(fid);
        

        
        %% do data sensitivity tests
     
        
        wtestsitedef.sites=[0 1e6 1e6];
        wtestsitedef.names={'GSL'};
        wtestsitedef.names2={'GSL'};
        wtestsitedef.firstage=min(oldest);

        wtestt=[1000 1800];
        Mdiff = [-1 1]/diff(wtestt);

        sitesub=(find(wdataset.siteid>9999));
        dosites=wdataset.siteid(sitesub);
        dositenames=wdataset.sitenames(sitesub);
        docoords=wdataset.sitecoords(sitesub,:);
        trainsub0=find((wdataset.limiting==0));
        
        clear sensf senssd sensnames senscoords;
        for qqq=0:length(dosites)
            for doexcl=0:1
                if qqq==0
                    if doexcl
                        trainsub=1;
                        wdataset2=wdataset; wdataset2.dY=ones(size(wdataset2.Y))*1e6; wdataset2.meantime=ones(size(wdataset2.Y))*1e5;
                        wdataset2.time1=wdataset2.meantime; wdataset2.time2=wdataset2.meantime;
                    else
                        wdataset2=wdataset; trainsub=trainsub0; Nincl(qqq+1)=length(trainsub);                   
                    end
                    sensnames{qqq+1} = 'All';
                    senscoords(qqq+1,:)=[NaN NaN];
                else
                    wdataset2=wdataset;
                    
                    if doexcl
                        trainsub=intersect(trainsub0,find(wdataset2.datid~=dosites(qqq)));
                    else                   
                        trainsub=intersect(trainsub0,find(wdataset2.datid==dosites(qqq))); Nincl(qqq+1)=length(trainsub);
                    end
                    sensnames{qqq+1}=dositenames{qqq};
                    senscoords(qqq+1,:)=docoords(qqq,:);
                end
                disp(sensnames{qqq+1});
                
                if length(trainsub)>0
                    
                     [wf,wsd,wV]=RegressHoloceneDataSets(wdataset2,wtestsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),wtestt,refyear,collinear);
                    wf2=Mdiff*wf; wsd2=sqrt(diag(Mdiff*wV*Mdiff'));
                else
                    wf2=NaN; wsd2=NaN;
                end
                
                
                sensf(qqq+1,doexcl+1) = wf2;
                senssd(qqq+1,doexcl+1) = wsd2;
                
            end
            
        end
        
        varnormalizer = senssd(1,2)^2;
        sensvarnormed = senssd.^2/varnormalizer;
        
        fid=fopen(['sitesensitivity' labl '.tsv'],'w');
        fprintf(fid,'Site sensitivity tests\n');
        fprintf(fid,'GSL rate, %0.0f to %0.0f\n\n',wtestt);
        fprintf(fid,'Site\tLat\tLong\tN\tInclusive f (mm/y)\t1s\tpercent var reduction\tpercent var reduction/N\tExclusive f (mm/y)\t1s\tpercent var increase\tpercent var increase/N\n');
        for qqq=1:size(sensf,1)
            fprintf(fid,sensnames{qqq});
            fprintf(fid,'\t%0.3f',senscoords(qqq,:));
            fprintf(fid,'\t%0.0f',Nincl(qqq));
            fprintf(fid,'\t%0.3f',[sensf(qqq,1) 1*senssd(qqq,1) 100*(1-sensvarnormed(qqq,1)) 100*(1-sensvarnormed(qqq,1))/Nincl(qqq) sensf(qqq,2) 1*senssd(qqq,2) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1)) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1))/Nincl(qqq)]);
            fprintf(fid,'\n');
        end        
        fclose(fid);
        
         
        %% do author sensitivity tests
        
        wtestsitedef.sites=[0 1e6 1e6];
        wtestsitedef.names={'GSL'};
        wtestsitedef.names2={'GSL'};
        wtestsitedef.firstage=min(oldest);

        wtestt=[1000 1800];
        Mdiff = [-1 1]/diff(wtestt);
        
        
        authors={'AUTHORTEAM','Gehrels','Gehrels-FirstAuthor','Horton','Engelhart','Long','Tornqvist','Kemp','Sivan'};
        studies={
                 {'EH12','Massachusetts','Spain','France','North Carolina','New Jersey','Florida','Connecticut'},
                 {'Isle of Wight','Maine','New Zealand','Tasmania','Nova Scotia','Scotland','Iceland'},
                 {'Maine','New Zealand','Tasmania','Nova Scotia','Iceland'},
                 {'EH12','Massachusetts-Wood Island','Spain','France','North Carolina','New Jersey','Florida'},
                 {'EH12'},
                 {'Scotland','Isle of Wight','Greenland'},
            {'Louisiana'},
            {'North Carolina','Florida','New Jersey','Massachusetts-Wood Island'},
            {'Israel'},
                };
        
        [s,si]=sort(authors);
        authors=authors(si); studies=studies(si);
        
        clear sitesubs dositess docoords;
        for qqq=1:length(authors)
            cursub=[];
            for rrr=1:length(studies{qqq})
                matches=strfind(wdataset.sitenames,studies{qqq}{rrr});
                for sss=1:length(matches)
                    if length(matches{sss})==1
                        cursub=union(cursub,sss);
                    end
                end
                
            end
            sitesubs{qqq}=cursub;
        end

        authors{end+1}='TIDEGAUGES';
        
        sitesubs{end+1}=find(ismember(wdataset.siteid,unique(wdataset.datid(find(wdataset.istg)))));
        
        for qqq=1:length(authors)
            dositess{qqq}=wdataset.siteid(sitesubs{qqq});
            docoords{qqq}=wdataset.sitecoords(sitesubs{qqq},:);
        end
        
        trainsub0=find((wdataset.limiting==0));
        
        clear sensf senssd sensnames senscoords;

        for qqq=0:length(dositess)
            for doexcl=0:1
                wdataset2=wdataset;
                if qqq==0
                    if doexcl
                        trainsub=1;
                        wdataset2=wdataset; wdataset2.dY=ones(size(wdataset2.Y))*1e6; wdataset2.meantime=ones(size(wdataset2.Y))*1e5;
                        wdataset2.time1=wdataset2.meantime; wdataset2.time2=wdataset2.meantime;
                    else
                        trainsub=trainsub0; Nincl(qqq+1)=length(trainsub);                   
                    end
                    sensnames{qqq+1} = 'ALL';
                else
                    
                    if doexcl
                        trainsub=intersect(trainsub0,find(~ismember(wdataset2.datid,dositess{qqq})));
                    else                   
                        trainsub=intersect(trainsub0,find(ismember(wdataset2.datid,dositess{qqq}))); Nincl(qqq+1)=length(trainsub);
                        trainsub=union(trainsub,find(wdataset2.istg));
                    end
                    sensnames{qqq+1}=authors{qqq};
                end
                disp(sensnames{qqq+1});
                
                % add in tide gauges
                
                if length(trainsub)>0
                    
                     [wf,wsd,wV]=RegressHoloceneDataSets(wdataset2,wtestsitedef,modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),wtestt,refyear,collinear);
                    wf2=Mdiff*wf; wsd2=sqrt(diag(Mdiff*wV*Mdiff'));
                else
                    wf2=NaN; wsd2=NaN;
                end
                
                
                sensf(qqq+1,doexcl+1) = wf2;
                senssd(qqq+1,doexcl+1) = wsd2;
                
            end
            
        end
        
        varnormalizer = senssd(1,2)^2;
        sensvarnormed = senssd.^2/varnormalizer;
        
        fid=fopen(['authorsensitivity' labl '.tsv'],'w');
        fprintf(fid,'Author sensitivity tests\n');
        fprintf(fid,'GSL rate, %0.0f to %0.0f\n\n',wtestt);
        fprintf(fid,'Author\tN\tInclusive f (mm/y)\t1s\tpercent var reduction\tpercent var reduction/N\tExclusive f (mm/y)\t1s\tpercent var increase\tpercent var increase/N\n');
        for qqq=1:size(sensf,1)
            fprintf(fid,sensnames{qqq});
            fprintf(fid,'\t%0.0f',Nincl(qqq));
            fprintf(fid,'\t%0.3f',[sensf(qqq,1) 1*senssd(qqq,1) 100*(1-sensvarnormed(qqq,1)) 100*(1-sensvarnormed(qqq,1))/Nincl(qqq) sensf(qqq,2) 1*senssd(qqq,2) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1)) 100*(sensvarnormed(qqq,2)-sensvarnormed(1,1))/Nincl(qqq)]);
            fprintf(fid,'\n');
        end        
        fclose(fid);
        
        
        
    
    save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','GISfpt','ICE5G','noiseMasks','testt','refyear');

end

% plot of relevant forcing proxies


clear proxyT proxyO;

proxyT.time1 = Chesapeake.yr;
sub=find((Chesapeake.yr>=1000).*(Chesapeake.yr<=1800));
Chesapeake.offset=mean(Chesapeake.T(sub));
proxyT.Y = (Chesapeake.T-Chesapeake.offset)/10;
proxyT.dY = ones(size(Chesapeake.T))*.01;
proxyT.datid = ones(size(Chesapeake.T));
proxyT.sitelen = length(Chesapeake.T);

proxyT.time1 = [proxyT.time1 ; round(GOM.yr)];
sub=find((GOM.yr>=1000).*(GOM.yr<=1800));
GOM.offset=mean(GOM.T(sub));
proxyT.Y = [proxyT.Y ; (GOM.T-GOM.offset)/5];
proxyT.dY = [proxyT.dY ; ones(size(GOM.T))*.01];
proxyT.datid = [proxyT.datid ; 2*ones(size(GOM.T))];
proxyT.sitelen = [proxyT.sitelen ; length(GOM.T)];

proxyT.time2 = proxyT.time1; proxyT.meantime=proxyT.time1;
proxyT.limiting = zeros(size(proxyT.Y));
proxyT.compactcorr = zeros(size(proxyT.Y));
proxyT.lat = proxyT.datid; proxyT.long = proxyT.datid;
proxyT.Ycv = sparse(diag(proxyT.dY).^2);
proxyT.istg = zeros(size(proxyT.Y));
proxyT.siteid = [1 2]';
proxyT.sitenames={'Chesapeake','GOM'};
proxyT.sitecoords=[1 1 ; 2 2];

%%

sub=find(RothJoos.yr>-1000);
proxyO.time1 = [RothJoos.yr(sub)];
proxyO.Y = [(RothJoos.TSI(sub) - 1365)*1000];
proxyO.dY = [ ones(size(RothJoos.TSI(sub)))*1];
proxyO.datid = [ones(size(RothJoos.yr(sub)))*1];
proxyO.sitelen = [ length(sub)];

sub=find((haug.yr>-1000).*(~isnan(haug.Ti)));
proxyO.time1 = [proxyO.time1 ; round(haug.yr(sub))];
proxyO.Y = [proxyO.Y ; haug.Ti(sub)*1000];
proxyO.dY = [proxyO.dY ; ones(size(haug.Ti(sub)))*1];
proxyO.datid = [proxyO.datid ; ones(size(haug.Ti(sub)))*2];
proxyO.sitelen = [proxyO.sitelen ; length(haug.Ti(sub))];

proxyO.time2 = proxyO.time1; proxyO.meantime=proxyO.time1;
proxyO.limiting = zeros(size(proxyO.Y));
proxyO.compactcorr = zeros(size(proxyO.Y));
proxyO.lat = proxyO.datid; proxyO.long = proxyO.datid;
proxyO.Ycv = sparse(diag(proxyO.dY).^2);
proxyO.istg = zeros(size(proxyO.Y));
proxyO.siteid = [1 2]';
proxyO.sitenames={'TSI','Cariacao'};
proxyO.sitecoords=[1 1 ; 2 2];


smoothwindow=101;

[proxyTs,proxyTthet] = GPSmoothTideGauges(proxyT,smoothwindow,2.0,20)
[proxyOs,proxyOthet] = GPSmoothTideGauges(proxyO,smoothwindow,2.0,20)
[proxyTs] = GPSmoothTideGauges(proxyT,smoothwindow,0,20,proxyTthet(1,:));

sub=find((proxyTs.datid==1).*(proxyTs.meantime>=(min(Chesapeake.yr)+smoothwindow/2)).*(proxyTs.meantime<=(max(Chesapeake.yr)-smoothwindow/2)));
ChesapeakeS.yr=proxyTs.meantime(sub);
ChesapeakeS.Tanom=proxyTs.Y(sub)*10;
ChespeakeS.T=ChesapeakeS.Tanom+Chesapeake.offset;

sub=find((proxyTs.datid==2).*(proxyTs.meantime>=(min(GOM.yr)+smoothwindow/2)).*(proxyTs.meantime<=(max(GOM.yr)-smoothwindow/2)));
GOMS.yr=proxyTs.meantime(sub);
GOMS.Tanom=proxyTs.Y(sub)*5;
GOMS.T=GOMS.Tanom+GOM.offset;

sub=find((proxyOs.datid==1).*(proxyOs.meantime>=(-1000+smoothwindow/2)).*(proxyOs.meantime<=(max(RothJoos.yr)-smoothwindow/2)));
RothJoosS.yr=proxyOs.meantime(sub);
RothJoosS.TSI=proxyOs.Y(sub)/1000+1365;

sub=find((proxyOs.datid==2).*(proxyOs.meantime>=(-1000+smoothwindow/2)).*(proxyOs.meantime<=(max(haug.yr)-smoothwindow/2)));
haugS.yr=proxyOs.meantime(sub);
haugS.Ti=proxyOs.Y(sub)/1000;

colrs={'b','g','r'};
clear ha;
figure;
hp=subplot(2,1,1);
ha(1)=plot(ChesapeakeS.yr,ChesapeakeS.Tanom,colrs{1},'linew',2); hold on
ha(2)=plot(GOMS.yr,GOMS.Tanom,colrs{2},'linew',2); hold on
ha(3)=plot(-1000,-1000,colrs{3},'linew',2);
set(hp,'box','off');
ylim([min([ChesapeakeS.Tanom ; GOMS.Tanom])-.5  max([ChesapeakeS.Tanom ; GOMS.Tanom])+.5]);
xlabel('Year');
hyl=ylabel('Sea Surface Temperature (C)'); %set(hyl,'Color',colrs{1});
legend(ha,sprintf('Chesapeake - %0.1f C',Chesapeake.offset),sprintf('Gulf of Mexico - %0.1f C',GOM.offset),'Florida Current (Sv)','Location','Southwest')

hp(2)=axes('Position',get(hp(1),'Position'));
ha(3)=plot(Lundtransport.yr,Lundtransport.Sv,colrs{3},'linew',2); hold on
plot(Lundtransport.yr,Lundtransport.Sv+Lundtransport.dSv,[colrs{3} '--']); hold on
plot(Lundtransport.yr,Lundtransport.Sv-Lundtransport.dSv,[colrs{3} '--']); hold on
set(hp(2),'Color','none','XAxisLoc','top','YAxisLoc','right','box','off','xtickl',{''});
ylim([27.5 32]);
hyl=ylabel({'Transport (Sv)'}); set(hyl,'Color',colrs{3});

set(hp,'xlim',[-500 2000]);
pdfwrite('forcingproxies1');

figure;

hp=subplot(2,1,1);
plot(RothJoosS.yr,RothJoosS.TSI,colrs{1},'linew',2); hold on
set(hp,'box','off','ylim',[8 17]);
sub=find(RothJoosS.yr>=-500);
ylim([min(RothJoosS.TSI(sub))-.05 max(RothJoosS.TSI(sub))+.05]);
hyl=ylabel('TSI (W/m^2)'); set(hyl,'Color',colrs{1});
xlabel('Year');

hp(2)=axes('Position',get(hp(1),'Position'));
plot(haugS.yr,haugS.Ti,colrs{2},'linew',1); hold on
set(hp(2),'Color','none','XAxisLoc','top','YAxisLoc','right','box','off','xtickl',{''});
sub=find(haugS.yr>=-500);
ylim([min(haugS.Ti(sub))-.02  max(haugS.Ti(sub))+.02]);
hyl=ylabel({'Cariaco Basin % Ti','(lower = dryer)'}); set(hyl,'Color',colrs{2});

set(hp,'xlim',[-500 2000]);
pdfwrite('forcingproxies2');

M=[ChesapeakeS.yr,ChesapeakeS.T];
save -ascii ChesapeakeS.asc M

M=[GOMS.yr,GOMS.T];
save -ascii GOM.asc M

M=[RothJoosS.yr,RothJoosS.TSI];
save -ascii RothJoosS.asc M

M=[haugS.yr,haugS.Ti];
save -ascii haugS.asc M
