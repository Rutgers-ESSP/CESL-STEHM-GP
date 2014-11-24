% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 14:58:44 EST 2014

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
        [fsl,Vsl]=SLRateMultisite(f2s{iii}(:,1),V2s{iii}(:,:,1),testsites(sitesub,:),testreg,testX(:,3),1000,1800);
        [fslA,sdslA,fslAd,sdslAd]=SLRateCompare(f2s{iii}(:,1),V2s{iii}(:,:,1),testsites(sitesub,:),testreg,testX(:,3),[0 900],[900 1800])

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
            wf=f2s{iii}(:,1); wV=V2s{iii}(:,:,1);
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
        
