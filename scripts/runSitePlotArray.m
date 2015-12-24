% Plot array elements for site plot.

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Dec 24 08:21:33 EST 2015


subsite=find((wdataset.siteid>10000));
maxdistfrom=0.1;
maxerror=1000;
wtestlocs=testlocs{iii};

% now array

% East River Marsh, Connecticut 20003
% Sand Point, North Carolina 140007
% Vioarholmi, Iceland 70001
% Loch Laxofrd, Scotland 190002
% Sissimut, Greenland 60003
% Caesarea, Israel 90001
% Christmas Island 10001
% Kariega Estuary, South Afirca 210002


letrs='abcdefgh';

sitesub = [20003 140007 70001 190002 60003 90001 10001 210002];
sitesubtitle={'East River Marsh, Connecticut','Sand Point, North Carolina','Vioarholmi, Iceland','Loch Laxford, Scotland','Sissimut, Greenland','Caesarea, Israel','Christmas Island, Kiribati',['Kariega Estuary, South Africa']};



clear testsitedefE;
testsitedefE.sites=[0 1e6 1e6];
testsitedefE.names={'GSL'};
testsitedefE.names2={'GSL'};
testsitedefE.firstage=min(oldest);
testsitedefE.oldest=min(oldest);
testsitedefE.youngest=2014;

for nnn=1:length(sitesub)
 sub=find(wdataset.datid==sitesub(nnn));
 testsitedefE.sites(end+1,:)=mean([wdataset.datid(sub) wdataset.lat(sub) wdataset.long(sub)],1);
 testsitedefE.names2={testsitedefE.names2{:}, sitesubtitle{nnn}};
 sublett=setdiff(1:length(sitesubtitle{nnn}),strfind(sitesubtitle{nnn},' '));
  sublett=setdiff(sublett,strfind(sitesubtitle{nnn},','));
 testsitedefE.names={testsitedefE.names{:}, sitesubtitle{nnn}(sublett)};
            testsitedefE.firstage = [testsitedefE.firstage min(oldest(:))];
            testsitedefE.oldest = [testsitedefE.oldest min(oldest(:))];
            testsitedefE.youngest = [testsitedefE.youngest 2014];
       testsitedefE.GIA = interp2(ICE5G.lat,ICE5G.long,ICE5G.gia,testsitedefE.sites(:,2),testsitedefE.sites(:,3),'linear');
    testsitedefE.GIA(find(testsitedefE.sites(:,2)>100))=0;
   
end

[fE,sdE,~,testlocsE]=RegressHoloceneDataSets(wdataset,testsitedefE,wmodelspec,thetTGG{jj},trainsub,noiseMasks,testt(subtimes),refyear,collinear);
fE=fE/10;
sdE=sdE/10;


for vvv=1:length(sitesub)
    kkk=find(testlocsE.sites(:,1)==sitesub(vvv));
    disp(sitesubtitle{vvv});
    
    clf;
    subplot(2,1,1);
    
    sub=find(testlocsE.reg==testlocsE.sites(kkk,1));
    if testlocsE.sites(kkk,2)<360
        distfrom=dDist(testlocsE.sites(kkk,2:3),[wdataset.lat wdataset.long]);
        subD=find(distfrom<maxdistfrom);
    else
        subD=find((wdataset.lat==testlocsE.sites(kkk,2)).*(wdataset.lat==testlocsE.sites(kkk,3)));
    end
    subD=intersect(subD,find(wdataset.dY<maxerror));
    
    clear plotdat;
    plotdat.x=testlocsE.X(sub,3);
    plotdat.y=fE(sub);
    plotdat.dy=sdE(sub)*2;
    
    PlotWithShadedErrors(plotdat,[0 0 0]);
    for uuu=subD(:)'
        plot([wdataset.time1(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)/10-2*wdataset.dY(uuu)/10*[1 1],'r'); hold on;
        plot([wdataset.time1(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)/10+2*wdataset.dY(uuu)/10*[1 1],'r');
        plot([wdataset.time1(uuu) wdataset.time1(uuu)],wdataset.Y0(uuu)/10+2*wdataset.dY(uuu)/10*[-1 1],'r');
        plot([wdataset.time2(uuu) wdataset.time2(uuu)],wdataset.Y0(uuu)/10+2*wdataset.dY(uuu)/10*[-1 1],'r');
    end
    plot(plotdat.x,plotdat.y,'k','linew',2);
    plot(plotdat.x,plotdat.y-plotdat.dy,'k--','linew',1);
    plot(plotdat.x,plotdat.y+plotdat.dy,'k--','linew',1);
    
    box on;
    set(gca,'fontsize',14);
    ht=title([letrs(vvv) ') ' sitesubtitle{vvv}]);
    set(ht,'fontsize',14);
    xlim([-500 2010]);
    ylabel('RSL (cm)');
    if vvv>=(length(sitesub)-1)
        xlabel('Year (CE)');
    else
        xlabel('   ');
    end
    pdfwrite(['siteplotarray-' letrs(vvv) '-' testlocsE.names{kkk}]);
    
end
