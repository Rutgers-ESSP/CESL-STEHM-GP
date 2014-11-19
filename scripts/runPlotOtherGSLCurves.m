% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 21:38:06 EST 2014

dat=importdata(fullfile(IFILES,'Grinsted2009_JonesA1B.txt'));
Grinsted2009_Jones.year=dat.data(:,1);
Grinsted2009_Jones.quantiles=[5 16 50 84 95];
Grinsted2009_Jones.y=dat.data(:,2:end)*1000;


dat=importdata(fullfile(IFILES,'Grinsted2009_MobergA1B.txt'));
Grinsted2009_Moberg.year=dat.data(:,1);
Grinsted2009_Moberg.quantiles=[5 16 50 84 95];
Grinsted2009_Moberg.y=dat.data(:,2:end)*1000;

% $$$ clf;
% $$$ refyear=2000;
% $$$ dodat=Grinsted2009_Moberg;
% $$$ sub=find(dodat.year==2000);
% $$$ doy = bsxfun(@minus,dodat.y,dodat.y(sub,:));
% $$$ plot(dodat.year,doy(:,3),'k','linew',2); hold on
% $$$ plot(dodat.year,doy(:,2),'k--'); hold on
% $$$ plot(dodat.year,doy(:,4),'k--'); hold on
% $$$ xlim([0 2010]);
% $$$ pdfwrite('Grinsted2009_Moberg');

% extract the detrended North Carolina sea-level curve

refyear=2000;
Mref = eye(size(testX,1));
for i=1:size(testsites,1)
    
    sub1=find(testreg==testsites(i,1));
    sub2=intersect(sub1,find(testX(:,3)==refyear));
    
    Mref(sub1,sub2)=Mref(sub1,sub2)-1;

end
Mref=sparse(Mref);

selsitenames={'NorthCarolina-SandPoint'};
sitesub=[];
for kk=1:length(selsitenames)
    q=find(strcmpi(selsitenames{kk},testsitedef.names));
    if length(q)>0
        sitesub=[sitesub q(1)];
    end
end
datsub=find(ismember(testreg,testsites(sitesub)));

colrs={'k'};

selmask=1;

wf=Mref*f2s{ii,jj}(:,selmask);
wV=Mref*V2s{ii,jj}(:,:,selmask)*Mref';
wsd=sqrt(diag(wV));


for dodetrend=[1]
    labl3='';
    if dodetrend
        [wf,wV,wsd]=DetrendSLReconstruction(wf,wV,testsites,testreg,testX(:,3),[0],1800,refyear);
        labl3='_detrended';
    end
    clf;
    [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,testsitedef.firstage(sitesub),testt(end),0,100,{'GSL'});
    set(hp,'xlim',[-500 2010]);
    
    labl2=[labl labl3];  

    %    delete(hp(2));
    %   pdfwrite(['NC_pseudoGSL' labl2]);
        
        NC_pseudoGSL=wf(datsub);
        NC_pseudoGSLsd=wsd(datsub);
        NC_yrs=testX(datsub,3);

end

%% do comparsion plot


clf; subplot(2,1,1); clear hp;
refyear=2000;
dodat=Grinsted2009_Moberg;
sub=find(dodat.year==2000);
sub1=find(dodat.year==0);
sub2=find(dodat.year==1800);

doy=dodat.y;
dorates=(doy(sub2,:)-doy(sub1,:))/1800;
doy=bsxfun(@minus,doy,bsxfun(@times,dorates(:,3),dodat.year-2000));
doy = bsxfun(@minus,doy,doy(sub,:));

hp(1)=plot(dodat.year,doy(:,3),'r','linew',2); hold on
plot(dodat.year,doy(:,2),'r--'); hold on
plot(dodat.year,doy(:,4),'r--'); hold on
xlim([0 2010]);

hold on;
hp(2)=plot(NC_yrs,NC_pseudoGSL,'g','linew',2); hold on;
plot(NC_yrs,NC_pseudoGSL+sqrt(NC_pseudoGSLsd.^2+50^2),'g--'); hold on;
plot(NC_yrs,NC_pseudoGSL-sqrt(NC_pseudoGSLsd.^2+50^2),'g--'); hold on;

datsub=find(testreg==0);
GSL=wf(datsub); GSLsd=wsd(datsub); GSL_yrs=testX(datsub,3);
hold on;
hp(3)=plot(GSL_yrs,GSL,'k','linew',3); hold on;
plot(GSL_yrs,GSL+GSLsd,'k--','linew',2); hold on;
plot(GSL_yrs,GSL-GSLsd,'k--','linew',2); hold on;

legend(hp,'Grinsted (Moberg)','North Carolina-derived','GSL','Location','Northwest');
ylabel('Detrended sea level (mm, \pm 1\sigma)'); ylim([-450 200]);
pdfwrite(['GSLcompare' labl]);
