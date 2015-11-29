% Compare GSL from statistical model to other sources.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 30 21:28:04 EST 2014

%% do comparsion plot

labl=labls{iii};
refyrs=[1600 1800];


for dodetrend=[0 1]
    datsub=find(testreg==0);

    dodat=Grinsted2009_Moberg;
    sub=find((dodat.year<=refyrs(2)).*(dodat.year>=refyrs(1)));
    sub1=find(dodat.year==0);
    sub2=find(dodat.year==1800);

    doy=dodat.y;
    dorates=(doy(sub2,:)-doy(sub1,:))/1800;
    if dodetrend
        doy=bsxfun(@minus,doy,bsxfun(@times,dorates(:,3),dodat.year-2000));
    end
    
    doy = bsxfun(@minus,doy,mean(doy(sub,3),1));

    clear plotdat;

    plotdat{1}.y=doy(:,3);
    plotdat{1}.dy=[doy(:,3)-doy(:,2) doy(:,4)-doy(:,3)];
    plotdat{1}.x=dodat.year;

    Mref=eye(length(NC_yrs));
    subt=find((NC_yrs<=refyrs(2)).*(NC_yrs>=refyrs(1)));
    Mref(:,subt)=Mref(:,subt)-1/length(subt);
    
    plotdat{2}.x=NC_yrs;
    plotdat{2}.y=Mref*NC_pseudoGSL;
    plotdat{2}.dy=sqrt(diag(Mref*NC_pseudoGSLV*Mref'));

    GSL=f2s{iii}(datsub,1); GSLV=V2s{iii}(datsub,datsub,1); GSL_yrs=testlocs{iii}.X(datsub,3);
    if dodetrend
        [GSL,GSLV,GSLsd]=DetrendSLReconstruction(GSL,GSLV,0,testlocs{iii}.reg(datsub),GSL_yrs,[0],1800,2000);
    end
    
    Mref=eye(length(GSL));
    subt=find((GSL_yrs<=refyrs(2)).*(GSL_yrs>=refyrs(1)));
    Mref(:,subt)=Mref(:,subt)-1/length(subt);

    plotdat{3}.y = Mref*GSL;
    plotdat{3}.dy = sqrt(diag(Mref*GSLV*Mref'));
    plotdat{3}.x = GSL_yrs;
    
    sub2000=find(GSL_yrs==2000);
    if ~dodetrend
        plotdat{4}.y=CW2011.y+plotdat{3}.y(sub2000)-CW2011.y(find(round(CW2011.year)==2000));
        plotdat{4}.dy=CW2011.dy;
        plotdat{4}.x=CW2011.year;
        
        plotdat{5}.x = Hay.meantime(:);
        plotdat{5}.y = Hay.Y(:)+plotdat{3}.y(sub2000)-Hay.Y(find(Hay.meantime==2000));
        plotdat{5}.dy = Hay.dY(:);
        
        rgb=[1 0 0 ; 0 1 0 ; 0 0 0 ; 1 1 0 ; 1 0 1];
        suborder=[1 2 4 5 3];
        plotlabls={'G09 hindcast','NC pseudo-GSL','GSL','CW2011','Hay2015'};
   else
        

        rgb=[1 0 0 ; 0 1 0 ; 0 0 0];
        suborder=1:3;
        plotlabls={'G09 hindcast','NC pseudo-GSL','GSL'};
    end
    

    if dodetrend
    clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[0 2010],[-450 200]);
        legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Detrended sea level (mm, \pm 1\sigma)'); box on;
        pdfwrite(['GSLcompare_detrended' labl]);
    else
    clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[0 2010],[-200 500]);
    legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Sea level (mm, \pm 1\sigma)'); box on;
       pdfwrite(['GSLcompare' labl]);
    
        clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[1800 2010],[-100 400]);
    legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Sea level (mm, \pm 1\sigma)'); box on;
       pdfwrite(['GSLcompare_1800' labl]);

    
    end
    
    
end
