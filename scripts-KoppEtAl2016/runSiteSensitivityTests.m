% Calculate sensitivity of results to different data subsets.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 01 11:41:15 EST 2016

%% do data sensitivity tests

outputGSLcurves=0;

% define intervals focused upon
firstyears=[0 700 1000 1400 1600 1800 1900];
lastyears= [700 1000 1400 1600 1800 1900 2000];
trainrange=[100 100 2000 2000];

wdat = datasets{ii};

% define subsets
sitesets0 = {'All','NWAtlantic','NEAtlantic','NAtlantic','SAmerica','AtlanticMediterranean'};
latbounds = [ -90 90 ; 0 90 ; 0 90 ; 0 90 ; -90 0 ; -90 90];
longbounds = [ -180 180 ; -125 -30 ; -30 30 ; -125 30 ; -125 0 ; -90 50];
sitelist = setdiff(wdat.siteid,0);
clear sitesetsub0;
for qqq=1:size(longbounds,1)
    sub=find((wdat.sitecoords(:,1)>=latbounds(qqq,1)).*(wdat.sitecoords(:,1)<=latbounds(qqq,2)));
    sub=intersect(sub,find(((wdat.sitecoords(:,2))>=longbounds(qqq,1)).*((wdat.sitecoords(:,2))<= longbounds(qqq,2))));
    sitesetsub0{qqq}=wdat.siteid(sub);
end

flattenersub=find((wdat.datid==0).*(wdat.dY>1e3));
clear sitesets sitesetsub;
for qqq=1:length(sitesets0)
    sitesets{qqq,1} = ['+' sitesets0{qqq} '+GSL'];
    sitesets{qqq,2} = ['+' sitesets0{qqq} '-GSL'];
    sitesets{qqq,3} = ['-' sitesets0{qqq} '+GSL'];
    sitesets{qqq,4} = ['-' sitesets0{qqq} '-GSL'];
    
    sitesetsub{qqq,2} = sitesetsub0{qqq};
    sitesetsub{qqq,1} = union(0,sitesetsub{qqq,2});
    sitesetsub{qqq,4} = setdiff(sitelist,sitesetsub0{qqq});
    sitesetsub{qqq,3} = union(0,sitesetsub{qqq,4});
    
end

ms=wmodelspec;
ms.thet0=thetTGG{jj}(1:end-1);
startcompact = thetTGG{jj}(end);

fid1=fopen(['sens_theta' labl '.tsv'],'w');
fid2=fopen(['sens_GSLrates' labl '.tsv'],'w');


for reoptimize=[0]

    for qqq=1:size(sitesetsub,1)
        for rrr=1:size(sitesetsub,2)
            
            disp(sitesets{qqq,rrr});
            
            datsub = find(ismember(wdat.datid,sitesetsub{qqq,rrr}));
            sitesub = find(ismember(wdat.siteid,sitesetsub{qqq,rrr}));
            
            if length(flattenersub)>0
                datsub=union(datsub,flattenersub);
                sitesub=union(sitesub,find(wdat.siteid==0));
            end
            
            wd=SubsetDataStructure(wdat,datsub,sitesub);
            
            dothet=thetTGG{jj};
            
            if length(datsub)>length(flattenersub)
                if reoptimize
                    [dothet]=OptimizeHoloceneCovariance(wd,ms,[3.0],trainfirsttime(end),trainrange(end),1e6,startcompact);  
                end
            else
                if length(flattenersub)==0
                    datsub=find(wdat.datid==0); datsub=datsub(end);
                    sitesub=find(wdat.siteid==0);
                    wd=SubsetDataStructure(wdat,datsub,sitesub);
                    wd.dY = 100000;
                    wd.Ycv = wd.dY^2;
                    wd.meantime=100000;
                    wd.time1=wd.meantime; wd.time2=wd.meantime;
                end
                
            end
            
            
            clear wdef;
            wdef.sites=[0 1e6 1e6];
            wdef.names={'GSL'};
            wdef.names2={'GSL'};
            wdef.firstage=min(firstyears);
            wdef.oldest=min(firstyears);
            wdef.youngest=2014;
            trainsub=find(wd.limiting==0);
	    if outputGSLcurves
                wtestt=testt;
	    else
                wtestt=min(firstyears):100:2000;
     end
     [wf,wsd,wV,wloc]=RegressHoloceneDataSets(wd,wdef,ms,dothet,trainsub,noiseMasks(1,:),wtestt,refyear,collinear);        

     [wfslope,wsdslope,wfslopediff,wsdslopediff,wdiffplus,wdiffless]=SLRateCompare(wf,wV,wloc.sites,wloc.reg,wloc.X(:,3),firstyears,lastyears);
     
     
     if (qqq==1) && (rrr==1)

         fprintf(fid1,'Trainset');
         fprintf(fid2,'Trainset');

         for pp=1:length(firstyears)
             fprintf(fid2,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
         end
         for pp=1:length(wdiffplus)
             fprintf(fid2,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
                           '%0.0f-%0.0f)\t2s'],[firstyears(wdiffplus(pp)) ...
                                 lastyears(wdiffplus(pp)) firstyears(wdiffless(pp)) lastyears(wdiffless(pp))]);
         end
         fprintf(fid2,'\n');
     end
     
     
     fprintf(fid1,sitesets{qqq,rrr});
     fprintf(fid2,sitesets{qqq,rrr});
     if reoptimize
         fprintf(fid1,'*');
         fprintf(fid2,'*');
     end
     
     for pp=1:length(firstyears)
         fprintf(fid2,'\t%0.2f',[wfslope(1,pp) 2*wsdslope(1,pp)]);
     end
     for pp=1:length(wdiffplus)
         fprintf(fid2,'\t%0.2f',[wfslopediff(1,pp) 2*wsdslopediff(1,pp)]);
     end
     
     fprintf(fid1,'\t%0.2f',dothet);
     
     fprintf(fid1,'\n');
     fprintf(fid2,'\n');
     
     sitesensfslope{iii,qqq,rrr,reoptimize+1}=wfslope(1,:);
     sitesenssdslope{iii,qqq,rrr,reoptimize+1}=wsdslope(1,:);
     sitesensfslopediff{iii,qqq,rrr,reoptimize+1}=wfslopediff(1,:);
     sitesenssdslopediff{iii,qqq,rrr,reoptimize+1}=wsdslopediff(1,:);
     sitesensthet{iii,qqq,rrr,reoptimize+1}=dothet;

     
     %%%%%%%%%%%%%%%
     if outputGSLcurves
         timesteps=100;

         sitesub=find(wloc.sites==0);
         datsub=find(wloc.reg==0);
         [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(wloc.X(datsub,3),wloc.reg(datsub),testsites(sitesub,1),wf,wV,colrs,wloc.X(datsub(1),3),testt(end),0,timesteps,{'GSL'});

         fid=fopen(['GSL_' labl '_sitesens_' sitesets{qqq,rrr} '.tsv'],'w');
         fprintf(fid,outtable);

         
         fclose(fid);
         
         
         % output GSL covariance

         fid=fopen(['GSL'  labl '_sitesens_' sitesets{qqq,rrr} '_cov.tsv'],'w');
         fprintf(fid,'mm^2');
         fprintf(fid,'\t%0.0f',testX(datsub,3));
         fprintf(fid,'\n');
         for ppp=1:length(datsub)
             fprintf(fid,'%0.0f',testX(datsub(ppp),3));
             fprintf(fid,'\t%0.8e',wV(datsub(ppp),datsub));
             fprintf(fid,'\n');
         end

         fclose(fid);
         
         
     end
     
     
     
        end
    end
end

fclose(fid1);
fclose(fid2);

sitesensfirstyears=firstyears;
sitesenslastyears=lastyears;
