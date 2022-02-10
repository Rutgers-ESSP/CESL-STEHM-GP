% Output table and plots of GSL.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, 2018-06-15 19:17:53 -0400

% operation matrices for zeroing to specific years

% operation matrices for zeroing to specific years

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
datsub=find(ismember(testreg,testsites(sitesub,1)));
colrs={'k'};

selmask=1;

%wf=Mref(datsub,datsub)*f2s{iii}(datsub,selmask);
%wV=Mref(datsub,datsub)*V2s{iii}(datsub,datsub,selmask)*Mref(datsub,datsub)';
wf=f2s{iii}(datsub,selmask);
wV=V2s{iii}(datsub,datsub,selmask);

wsd=sqrt(diag(wV));

wf1900=Mref2(datsub,datsub)*f2s{iii}(datsub,selmask);
wV1900=Mref2(datsub,datsub)*V2s{iii}(datsub,datsub,selmask)*Mref2(datsub,datsub)';
wsd1900=sqrt(diag(wV1900));


for dodetrend=[0]
    labl3='';
    
% $$$     if dodetrend
% $$$         %        [wf,wV,wsd]=DetrendSLReconstruction(wf,wV,testsites(sitesub),testreg(datsub),testX(datsub,3),[0],1800,refyear);
% $$$         % labl3='_detrended';
% $$$         labl3='_nonegcov';
% $$$         wV=wV.*(wV>0);
% $$$     end
    labl2=[labls{iii} labl3];
    
    if ~dodetrend
        
        
        clear plotdat;
        plotdat{1}.y = wf;
        plotdat{1}.dy = sqrt(diag(wV));
        plotdat{1}.x = testX(datsub,3); 
        rgb=[0 0 0];
        
        clf;
        subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat,rgb,[],[],[],[-1000 2010]);
        if dodetrend
            ylabel('Detrended GSL (mm, \pm 1\sigma)'); 
        else
            ylabel('GSL (mm, \pm 1\sigma)');
        end
        

        pdfwrite(['GSL' labl2]);
        
    end
    

    for timesteps=[100 200 400 40 20]
        [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf,wV,colrs,testsitedef.firstage(sitesub),testt(end),0,timesteps,{'GSL'});
        set(hp,'xlim',[-1000 2010]);
        delete(hp(2));
        
        clf;
        
        if (~dodetrend)
            clear plotdat;
            plotdat{1}.y = dGSL;
            plotdat{1}.dy = dGSLsd;
            plotdat{1}.x = difftimes; 
            rgb=[0 0 0];
            
            subplot(2,1,1);
            plot([-1000 2010],[0 0],'k-'); hold on;
            [hl,hk]=PlotWithShadedErrors(plotdat,rgb,[],[],[],[-1000 2010]);
            ylabel([num2str(timesteps) '-year avg. GSL rate (mm/y, \pm 1\sigma)']);
            pdfwrite(['GSLrate_' num2str(timesteps) 'y_' labl2]);
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
        wquants=[.67 .9 .95 .96 .97 .98 .99 .995 .999];
        if length(subyrs)>0
            fprintf(fid,'\n\nLast year faster');
            fprintf(fid,'\t%0.3f',wquants);
            fprintf(fid,'\n');
            
            clear lastyrlarger;
            for mm=1:length(subyrs)
                curyr=difftimes(suboverallyrs(subyrs(mm)));
                compsamps=samps(:,suboverallyrs(subyrs(mm)));
                sub=intersect(suboverallyrs,find(difftimes<curyr));
                largeryrs=bsxfun(@gt,samps(:,sub),compsamps);
                largeryrs2=bsxfun(@times,largeryrs,difftimes(sub)') -1e6*(~largeryrs);
                largeryrs2(find(largeryrs2<-1e5))=min(difftimes)-1;
                lastyrlarger(:,mm)=max(largeryrs2,[],2);
                fprintf(fid,'\n%0.0f',curyr);
                fprintf(fid,'\t%0.0f',quantile(lastyrlarger(:,mm),wquants));
            end
        end
        fclose(fid);
        
        
    end    
    % output GSL covariance


    if ~dodetrend
        
        fid=fopen(['GSL'  labl2 '_cov.tsv'],'w');
        fprintf(fid,'mm^2');
        fprintf(fid,'\t%0.0f',testX(datsub,3));
        fprintf(fid,'\n');
        for ppp=1:length(datsub)
            fprintf(fid,'%0.0f',testX(datsub(ppp),3));
            fprintf(fid,'\t%0.8e',wV(datsub(ppp),datsub));
            fprintf(fid,'\n');
        end

        fclose(fid);
        
        clf;
        wcorr=wV./bsxfun(@times,sqrt(diag(wV)),sqrt(diag(wV))');
        imagesc(testX(datsub,3),testX(datsub,3),wcorr(datsub,datsub));
        colorbar;
        pdfwrite(['GSL'  labl2 '_corr']);

    end
    
% $$$     [wVNL,wVconst,wVrate,wt0,wsigadj]=PartitionCovarianceLNL(wV,testX(datsub,3));
% $$$     fid=fopen(['GSL'  labl2 '_covNL.tsv'],'w');
% $$$     fprintf(fid,'Vconst = %0.8e\n',wVconst);
% $$$     fprintf(fid,'Vrate = %0.8e\n',wVrate);
% $$$     fprintf(fid,'t0 = %0.8e\n',wt0);
% $$$     fprintf(fid,'adj = %0.8e\n',wsigadj);
% $$$     fprintf(fid,'\n\n',wVconst);
% $$$ 
% $$$     fprintf(fid,'mm^2');
% $$$     fprintf(fid,'\t%0.0f',testX(datsub,3));
% $$$     fprintf(fid,'\n');
% $$$     for ppp=1:length(datsub)
% $$$         fprintf(fid,'%0.0f',testX(datsub(ppp),3));
% $$$         fprintf(fid,'\t%0.8e',wVNL(datsub(ppp),datsub));
% $$$         fprintf(fid,'\n');
% $$$     end
% $$$ 
% $$$     fclose(fid);

    fid=fopen(['dGSL_' num2str(timesteps) 'y' labl2 '_cov.tsv'],'w');
    fprintf(fid,'(mm/y)^2');
    fprintf(fid,'\t%0.0f',difftimes);
    fprintf(fid,'\n');
    for ppp=1:length(difftimes)
        fprintf(fid,'%0.0f',difftimes(ppp));
        fprintf(fid,'\t%0.8e',dGSLV(ppp,:));
        fprintf(fid,'\n');
    end
    fclose(fid);   
end