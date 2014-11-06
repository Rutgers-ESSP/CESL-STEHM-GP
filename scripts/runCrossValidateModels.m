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