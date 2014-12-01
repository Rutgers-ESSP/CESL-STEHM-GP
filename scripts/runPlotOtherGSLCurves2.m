% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 30 21:20:32 EST 2014

%% do comparsion plot

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

    GSL=f2s{1}(datsub,1); GSLV=V2s{1}(datsub,datsub,1); GSL_yrs=testlocs{1}.X(datsub,3);
    if dodetrend
        [GSL,GSLV,GSLsd]=DetrendSLReconstruction(GSL,GSLV,0,testlocs{1}.reg(datsub),GSL_yrs,[0],1800,2000);
    end
    
    Mref=eye(length(GSL));
    subt=find((GSL_yrs<=refyrs(2)).*(GSL_yrs>=refyrs(1)));
    Mref(:,subt)=Mref(:,subt)-1/length(subt);

    plotdat{3}.y = Mref*GSL;
    plotdat{3}.dy = sqrt(diag(Mref*GSLV*Mref'));
    plotdat{3}.x = GSL_yrs;

    sub2000=find(GSL_yrs==2000);

    
    GSL2=f2s{2}(datsub,1); GSLV2=V2s{2}(datsub,datsub,1);
    if dodetrend
        [GSL2,GSLV2,GSLsd2]=DetrendSLReconstruction(GSL2,GSLV2,0,testlocs{2}.reg(datsub),GSL_yrs,[0],1800,2000);
    end
    
    plotdat{4}.y = Mref*GSL2;
    %plotdat{4}.y = plotdat{4}.y-plotdat{4}.y(sub2000)+plotdat{3}.y(sub2000);
    plotdat{4}.dy = sqrt(diag(Mref*GSLV2*Mref'));
    plotdat{4}.x = GSL_yrs;    
    
    if ~dodetrend
        plotdat{5}.y=CW2011.y+plotdat{3}.y(sub2000)-CW2011.y(find(round(CW2011.year)==2000));
        plotdat{5}.dy=CW2011.dy;
        plotdat{5}.x=CW2011.year;
        
        plotdat{6}.x = Hay.meantime(:);
        plotdat{6}.y = Hay.Y(:)+plotdat{3}.y(sub2000)-Hay.Y(find(Hay.meantime==2000));
        plotdat{6}.dy = Hay.dY(:);
        
        rgb=[1 0 0 ; 0 1 0 ; 0 0 0 ; 1 .5 0 ; 1 1 0 ; 1 0 1];
        suborder=[1 2 5 6 4 3];
        plotlabls={'G09','NC','GSL_{ML}','GSL_{Gr}','CW2011','Hay2015'};
   else
        

        rgb=[1 0 0 ; 0 1 0 ; 0 0 0; 1 .5 0];
        suborder=[1 2 4 3];
        plotlabls={'G09 hindcast','NC pseudo-GSL','GSL_{ML}','GSL_{Gr}'};
    end
    

    if dodetrend
    clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[-1000 2010],[-450 200]);
        legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Detrended sea level (mm, \pm 1\sigma)'); box on;
        pdfwrite(['GSLcompare_detrended']);
    else
    clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[-1000 2010],[-200 500]);
    legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Sea level (mm, \pm 1\sigma)'); box on;
       pdfwrite(['GSLcompare']);
    
        clf;
    subplot(2,1,1);
        [hl,hk]=PlotWithShadedErrors(plotdat(suborder),rgb(suborder,:),[],[],[],[1800 2010],[-100 400]);
    legend(hl,plotlabls(suborder),'Location','Northwest');
        ylabel('Sea level (mm, \pm 1\sigma)'); box on;
       pdfwrite(['GSLcompare_1800']);

    
    end
    
    
end
