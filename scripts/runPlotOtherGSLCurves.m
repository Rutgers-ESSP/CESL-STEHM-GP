% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 21:38:06 EST 2014

dat=importdata(fullfile(IFILES,'Grinsted2009_JonesA1B.txt'));
Grinsted2009_Jones.year=dat.data(:,1);
Grinsted2009_Jones.quantiles=[5 16 50 84 95];
Grinsted2009_Jones.y=dat.data(:,2:end);


dat=importdata(fullfile(IFILES,'Grinsted2009_MobergA1B.txt'));
Grinsted2009_Moberg.year=dat.data(:,1);
Grinsted2009_Moberg.quantiles=[5 16 50 84 95];
Grinsted2009_Moberg.y=dat.data(:,2:end);

clf;
refyear=2000;
dodat=Grinsted2009_Moberg;
sub=find(dodat.year==2000);
doy = bsxfun(@minus,dodat.y,dodat.y(sub,:));
plot(dodat.year,doy(:,3),'k','linew',2); hold on
plot(dodat.year,doy(:,2),'k--'); hold on
plot(dodat.year,doy(:,4),'k--'); hold on
xlim([0 2010]);
pdfwrite('Grinsted2009_Moberg');

% extract the detrended North Carolina sea-level curve

refyear=2000;
Mref = eye(size(testX,1));
for i=1:size(testsites,1)
    
    sub1=find(testreg==testsites(i,1));
    sub2=intersect(sub1,find(testX(:,3)==refyear));
    
    Mref(sub1,sub2)=Mref(sub1,sub2)-1;

end
Mref=sparse(Mref);

selsitenames={'NorthCarolina-TumpPoint'};
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

    delete(hp(2));
    if timesteps==100
        pdfwrite(['NC_pseudoGSL' labl2]);
    end

end