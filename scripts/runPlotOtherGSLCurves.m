% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 20 19:28:20 EST 2014

%% do comparsion plot

datsub=find(testreg==0);

refyear=2000;
dodat=Grinsted2009_Moberg;
sub=find(dodat.year==2000);
sub1=find(dodat.year==0);
sub2=find(dodat.year==1800);

doy=dodat.y;
dorates=(doy(sub2,:)-doy(sub1,:))/1800;
doy=bsxfun(@minus,doy,bsxfun(@times,dorates(:,3),dodat.year-2000));
doy = bsxfun(@minus,doy,doy(sub,:));

clear plotdat;

plotdat{1}.y=doy(:,3);
plotdat{1}.dy=[doy(:,3)-doy(:,2) doy(:,4)-doy(:,3)];
plotdat{1}.x=dodat.year;

plotdat{2}.x=NC_yrs;
plotdat{2}.y=NC_pseudoGSL;
plotdat{2}.dy=NC_pseudoGSLsd;

GSL=f2s{iii}(datsub,1); GSLV=V2s{iii}(datsub,datsub,1); GSL_yrs=testlocs{iii}.X(datsub,3);
[GSL,GSLV,GSLsd]=DetrendSLReconstruction(GSL,GSLV,0,testlocs{iii}.reg(datsub),GSL_yrs,[0],1800,2000);
Mref=eye(length(GSL));
subt=find(GSL_yrs==2000);
Mref(:,subt)=Mref(:,subt)-1;

plotdat{3}.y = Mref*GSL;
plotdat{3}.dy = sqrt(diag(Mref*GSLV*Mref'));
plotdat{3}.x = GSL_yrs; 

rgb=[1 0 0 ; 0 1 0 ; 0 0 0];
suborder=1:3;
plotlabls={'G09 hindcast','NC pseudo-GSL','GSL'};

clf;
subplot(2,1,1);
[hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[0 2010],[-450 200]);

%xlim([0 2010]); ylim([-450 200]);
legend(hl,plotlabls(suborder),'Location','Northwest');
ylabel('Detrended sea level (mm, \pm 1\sigma)'); 
box on;
pdfwrite(['GSLcompare' labl]);
