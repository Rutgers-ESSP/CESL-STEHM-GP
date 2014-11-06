% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 14:54:30 EST 2014


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
            [wf,wV,wsd]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[0],1800,refyear);
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
   
    