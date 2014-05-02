% Master script for Common Era proxy + tide gauge sea-level analysis
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Apr 30 13:47:30 EDT 2014

pd=pwd;
addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pd '/MFILES']);
IFILES=[pd '/IFILES'];
addpath(pd)
savefile='~/tmp/CESL';

WORKDIR='140430b';
if ~exist(WORKDIR,'dir')
    mkdir(WORKDIR);
end
cd(WORKDIR);


%
firsttime=-1000;

runSetupHoloceneCovariance;
runImportHoloceneDataSets;

addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pd '/MFILES']);
addpath(pd)

save(savefile,'datasets','modelspec');
addpath(pwd)
trainsets=[1 2 ]; % train w/TG+GSL+PX or w/TG+PX
trainspecs=[1 1];
trainlabels={'default','trainedwoGSL'};

clear thetTGG trainsubsubset;
for ii=1:length(trainspecs)
    [thetTGG{ii},trainsubsubset{ii},logp(ii)]= ...
        OptimizeHoloceneCovariance(datasets{trainsets(ii)}, ...
                                   modelspec(trainspecs(ii)),[2.4 2.0 3.0]);
    
    fid=fopen('thetTGG.tsv','w');
    for iii=1:length(thetTGG)
        fprintf(fid,[datasets{trainsets(iii)}.label '\t' ...
                     modelspec(trainspecs(iii)).label '\t']);
        fprintf(fid,['(%0.2f)\t'],logp(iii));
        fprintf(fid,'%0.3f\t',thetTGG{iii});
        fprintf(fid,'\n');
    end
    fclose(fid);
end

save thetTGG thetTGG trainsubsubset

%%%%%%%%


%%%%

fid=fopen('TGandProxyData.tsv','w');
fprintf(fid,'ID\t1ime1\ttime2\tlimiting\tY-GIAproj\tY\tdY\tcompactcorr\tistg\tlat\tlong\n');
ii=1;
for i=1:size(datasets{ii}.datid,1)
    compactcorr=full(datasets{ii}.compactcorr);
    fprintf(fid,'%d\t',datasets{ii}.datid(i));
    fprintf(fid,'%0.1f\t',datasets{ii}.time1(i));
    fprintf(fid,'%0.1f\t',datasets{ii}.time2(i));
    fprintf(fid,'%d\t',datasets{ii}.limiting(i));
    fprintf(fid,'%0.2f\t',datasets{ii}.Y(i));
    fprintf(fid,'%0.2f\t',datasets{ii}.Y0(i));
    fprintf(fid,'%0.2f\t',datasets{ii}.dY(i));
    fprintf(fid,'%0.2f\t',compactcorr(i));
    fprintf(fid,'%d\t',datasets{ii}.istg(i));
    fprintf(fid,'%0.2f\t',datasets{ii}.lat(i));
    fprintf(fid,'%0.2f\n',datasets{ii}.long(i));
end
fclose(fid)

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

regresssets=[2 1 2 3 4 5 6 7]; % regress w/TG+PX
regressparams=[1 1 2 1 1 1 1 1];
clear regresslabels;
for i=1:length(regresssets)
    regresslabels{i} = [datasets{regresssets(i)}.label '_' trainlabels{regressparams(i)}];
end

for iii=1:length(regresssets)
    ii=regresssets(iii);
    jj=regressparams(iii);

    noiseMasks = ones(6,length(thetTGG{trainspecs(jj)}));
    noiseMasks(1,[modelspec(trainspecs(jj)).subampnoise])=0; % without white noise
    noiseMasks(2,[modelspec(trainspecs(jj)).subamplinear modelspec(trainspecs(jj)).subampoffset]  )=0; %without linear
    noiseMasks(3,[modelspec(trainspecs(jj)).subampglobal modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local
    noiseMasks(4,[modelspec(trainspecs(jj)).subamplinear modelspec(trainspecs(jj)).subampglobal  modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local non-linear
    noiseMasks(5,[modelspec(trainspecs(jj)).subampglobal modelspec(trainspecs(jj)).subampnonlinear modelspec(trainspecs(jj)).subampoffset])=0; %only regional and local linear
    noiseMasks(6,setdiff(modelspec(trainspecs(jj)).subamp,modelspec(trainspecs(jj)).subampGIS))=0; %Greenland only
    noiseMasklabels={'denoised','nonlin','regloc','reglocnonlin','regloclin','gis'};

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
    [fp(subps),sdp(subps),~,testlocp]=RegressHoloceneDataSets(wdataset,testsitedefp,                                            modelspec(trainspecs(jj)),thetTGG{jj},trainsub,noiseMasks(1,:),testtsp,refyear,collinear,passderivs,invcv);

    fid=fopen(['TGandProxyData' labl '.tsv'],'w');
    fprintf(fid,['Site\tID\t1ime1\ttime2\tlimiting\tGIAproj\tY-GIAproj\tY\' ...
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
    
    makeplots_sldecomp(wdataset,f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},labl);
    
    testreg=testlocs{ii,jj}.reg;
    testsites=testlocs{ii,jj}.sites;
    testX=testlocs{ii,jj}.X;
    testnames2=testlocs{ii,jj}.names2;

    firstyears=[ 0 1000  200  -300 700 1700 1800 1900 860  1000 1080 1460 1560 1800 1840];
    lastyears=[1800 1800 1000 700 1700 1800 1900 2000 1560 1400 1360 1880 1800 1880 1900];
    
    clear fslopelin fslopeavg sdslopelin sdslopeavg fslopeavgdiff sdslopeavgdiff;

    for kk=1:size(testsites,1)
        
        %%%
        sub=find((testreg==testsites(kk,1)));
        M=zeros(length(firstyears),length(sub));
        for pp=1:length(firstyears)

            sub1=find((testX(sub,3)==firstyears(pp)));
            sub2=find((testX(sub,3)==lastyears(pp)));
            if (length(sub1)==1)&&(length(sub2)==1)
                M(pp,sub1(1))=-1; M(pp,sub2(1))=1;
                M(pp,:)=M(pp,:)/(testX(sub2(1),3)-testX(sub1(1),3));
            end
        end
        
        clear Mdiff diffplus diffless;
        counter=1;
        for pp=1:min(7,length(firstyears))
            for qq=(pp+1):length(firstyears)
                Mdiff(counter,pp)=-1;
                Mdiff(counter,qq)=1;
                diffplus(counter)=qq; diffless(counter)=pp;
                counter=counter+1;
            end
        end
        
        fslopelin(kk,ii,jj) = M(end,:)*f2s{ii,jj}(sub,5);
        Vslope = M(end,:)*V2s{ii,jj}(sub,sub,5)*M(end,:)';
        sdslopelin(kk,ii,jj)=sqrt(diag(Vslope));
        
        
        fslope=M*f2s{ii,jj}(sub,1);
        Vslope=M*V2s{ii,jj}(sub,sub,1)*M';
        sdslope=sqrt(diag(Vslope));

        notgood=find(sum(abs(M),2)==0);
        modifier=0*sdslope; modifier(notgood)=1e6;
        sdslope=sdslope+modifier; Vslope=Vslope+diag(modifier.^2);

        fslopediff = Mdiff*fslope;
        Vslopediff = Mdiff*Vslope*Mdiff';
        sdslopediff = sqrt(diag(Vslopediff));
        
        fslope(find(sdslope>1e4))=NaN;
        sdslope(find(sdslope>1e4))=NaN;
        fslopediff(find(sdslopediff>1e4))=NaN;
        sdslopediff(find(sdslopediff>1e4))=NaN;
        

        fslopeavg(kk,ii,jj,:) = reshape(fslope,1,1,1,[]);
        sdslopeavg(kk,ii,jj,:)= reshape(sdslope,1,1,1,[]);

        fslopeavgdiff(kk,ii,jj,:) = reshape(fslopediff,1,1,1,[]);
        sdslopeavgdiff(kk,ii,jj,:)= reshape(sdslopediff,1,1,1,[]);

 
    end
    
    
    fid=fopen(['linrates' labl '.tsv'],'w');
    fprintf(fid,['Regional+Local Linear Rates (mm/y), ' labl '\n']);
    fprintf(fid,'Site\tICE5G VM2-90\tRate (linear)\t2s');
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
        fprintf(fid,'\t%0.2f',testsitedef.GIA(kk));
        fprintf(fid,'\t%0.2f',[fslopelin(kk,ii,jj) 2*sdslopelin(kk,ii,jj)]);
        for pp=1:length(firstyears)
            fprintf(fid,'\t%0.2f',[fslopeavg(kk,ii,jj,pp) 2*sdslopeavg(kk,ii,jj,pp)]);
        end
        for pp=1:length(diffplus)
            fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,ii,jj,pp) 2*sdslopeavgdiff(kk,ii,jj,pp)]);
        end
        fprintf(fid,'\n');
    end
    fclose(fid);

    
    fslope=fslopeavg(:,ii,jj,2); fslope=fslope(:);
    sub=find(testsites(:,2)<1e3);

    figure;
    plotcont; hold on;

    % plot tide gauge locations
    scatter(mod(TGNOCW.sitecoords(:,2),360),TGNOCW.sitecoords(:,1),10,[.5 .5 .5],'filled'); hold on;

    % plot observation sites
    scatter(mod(testsites(sub,3),360),testsites(sub,2),30,fslope(sub),'filled'); hold on;

    axis tight;
    colorbar;
    box on;
    %title('Sea level  rates, 1000-1800 CE (mm/y)');
    pdfwrite(['map_linrates_global' labl]);

    figure;
    usstates = shaperead('usastatehi', 'UseGeoCoords', true);
    plot((mod([usstates.Lon],360)),[usstates.Lat],'k','linew',.5); hold on;
    if exist([IFILES '/province/PROVINCE.SHP']); canada = shaperead([IFILES '/province/PROVINCE.SHP'],'UseGeoCoords',true); else; canada=[]; end
    if length(canada)>0 ; plot((mod([canada.Lon],360)),[canada.Lat],'k','linew',.5); end
    xl=xlim;
    yl=ylim;

    % plot tide gauge locations
    scatter(mod(TGNOCW.sitecoords(:,2),360),TGNOCW.sitecoords(:,1),30,[.3 .3 .3],'filled'); hold on;

    % plot observation sites
    scatter(mod(testsites(sub,3),360),testsites(sub,2),100,fslope(sub),'filled'); hold on;

    axis tight;
    colorbar;
    box on;
    %title('Sea level  rates, 1000-1800 CE (mm/y)');

    xlim([260 305]);
    ylim([27 50]);

    pdfwrite(['map_linrates_NA' labl]);


    %%% 

    refyear=2000;
    Mref = eye(size(testX,1));
    for i=1:size(testsites,1)
        
        sub1=find(testreg==testsites(i,1));
        sub2=intersect(sub1,find(testX(:,3)==refyear));
        
        Mref(sub1,sub2)=Mref(sub1,sub2)-1;

    end
    Mref=sparse(Mref);
    selsitenames={'Florida-Nassau','NorthCarolina-TumpPoint','NorthCarolina-SandPoint','NewJersey-LeedsPoint','NovaScotia-Chezzetcook'};
    shortnames = {'FL','NC-TP','NC-SP','NJ','NS'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end
    colrs={'r','c','b','m','g'};

    datsub=find(ismember(testreg,testsites(sitesub,1)));

    % figure of sites with data, full curves
    % figure of sites with data, without linear

    for selmask=[1 2]
        
        clf;
        
        wf=Mref*f2s{ii,jj}(:,selmask);
        wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
        wsd=sqrt(diag(wV));
        [hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wsd(datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,100,testnames2(sitesub));
        
        if selmask==1
            legend(hl,shortnames,'Location','Southwest');
        else
            legend(hl,shortnames,'Location','Southwest');
        end
        set(hp,'xlim',[-1000 2010]);
        %	set(hp(2),'ylim',[-2 5]);

        pdfwrite(['paleorsl_' noiseMasklabels{selmask} labl]);

        fid=fopen(['paleorsl_' noiseMasklabels{selmask} labl '.tsv'],'w');
        fprintf(fid,outtable);
        fclose(fid);
    end


    % figure of regional averages, without linear or GSL

    sitesub=sitesub(2:end);
    colrs=colrs(2:end);
    datsub=find(ismember(testreg,testsites(sitesub,1)));


    for selmask=[3 4]

        wf=Mref*f2s{ii,jj}(:,selmask);
        wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
        wsd=sqrt(diag(wV));
        
        clf;
        [hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wsd(datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,160,testnames2(sitesub));

        legend(hl,testnames2{sitesub},'Location','Southwest');

        set(hp,'xlim',[-1000 2010]);

        pdfwrite(['paleorsl_' noiseMasklabels{selmask} labl]);

        fid=fopen(['paleorsl_' noiseMasklabels{selmask} labl '.tsv'],'w');
        fprintf(fid,outtable);
        fclose(fid);
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
    colrs={'k','r','c','b','g','m','y'};

    selmask=1;

    wf=Mref*f2s{ii,jj}(:,selmask);
    wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
    wsd=sqrt(diag(wV));

    datsub=find(ismember(testreg,testsites(sitesub,1)));

    dorateshift=0;
    if length(intersect(modelspec(trainspecs(jj)).subamplinear, ...
                        modelspec(trainspecs(jj)).subampglobal))==0
        dorateshift=[0 0.15 0.3];
    elseif sum(abs(intersect(modelspec(trainspecs(jj)).subamplinear, ...
                             modelspec(trainspecs(jj)).subampglobal))) == 0
        dorateshift=[0 0.15 0.3];
    end
    
    
    for rateshift=dorateshift
        for timesteps=[100 1000 160 40 20 1800]

            clf;
            [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub)+rateshift*(testX(datsub,3)-refyear),wsd(datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,timesteps,{'GSL'});
            set(hp,'xlim',[-1000 2010]);
            
            labl2=labl;
            if abs(rateshift)~=0
                labl2=[labl2 '_shift' sprintf('%0.0f',rateshift*100)];
            end
            

            pdfwrite(['GSL_' num2str(timesteps) 'y' labl2]);



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
    end
    
    
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
        [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wsd(datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,timesteps,{'GIS'});
        set(hp,'xlim',[-1000 2010]);

        pdfwrite(['GIS_' num2str(timesteps) 'y' labl]);


        fid=fopen(['GIS_' num2str(timesteps) 'y' labl '.tsv'],'w');
        fprintf(fid,outtable);
        
    end

    %%%%

    % figure of sea-level differences

    selsitenames={'Florida-Nassau','NorthCarolina-SandPoint','NewJersey-LeedsPoint','Massachusetts','NovaScotia-Chezzetcook'};
    shortselnames={'FL','NC','NJ','MA','NS'};
    sitesub=[];
    for kk=1:length(selsitenames)
        q=find(strcmpi(selsitenames{kk},testsitedef.names));
        if length(q)>0
            sitesub=[sitesub q(1)];
        end
    end

    selmask=2;

    clear pairsets;
    pairsets{1}=[2 1; 2 3; 2 4; 2 5 ; 3 5];
    pairsets{2}=[2 3; 3 5 ; 2 5];
    for mmm=1:length(pairsets)
        gradf=[]; gradsd=[]; gradt=[]; gradpair=[]; gradV=[]; ...
              gradstarttimes=[]; clear pairnames;
        diffpairs=pairsets{mmm};
        labl2=['_' num2str(mmm)];
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

            gradf = [gradf ; Mgrad*f2s{ii,jj}(:,selmask)];
            gradV(length(gradV) + [1:length(u)],length(gradV) + [1:length(u)]) = Mgrad*V2s{ii,jj}(:,:,selmask)*Mgrad';
            gradt = [gradt ; u];
            gradpair = [gradpair ; ones(length(u),1)*i];
            gradstarttimes(i)=u(1);
            
        end 
        gradsd = sqrt(diag(gradV));

        for timesteps=[100 600 160 40 20]

            clf;
            [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(gradt,gradpair,1:length(diffpairs),gradf,gradV,colrs,gradstarttimes,2010,0,timesteps,pairnames);
            set(hp,'xlim',[-1000 2010]);
            if mmm>1
                set(hp,'xlim',[0 2010]);
            end
            legend(hl,pairnames,'Location','Southwest');
            pdfwrite(['RSLgrad_' num2str(timesteps) 'y' labl labl2]);

            if mmm==1
                fid=fopen(['RSLgrad_' num2str(timesteps) 'y' labl '.tsv'],'w');
                fprintf(fid,outtable);
                fclose(fid);
            end
            
        end
    end
    
    save(savefile,'datasets','modelspec','f2s','sd2s','V2s', ...
         'testlocs','logp','testsitedef','trainspecs','thetTGG','GISfpt','ICE5G','noiseMasks','testt','refyear');

end
