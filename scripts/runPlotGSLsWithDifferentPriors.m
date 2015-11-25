% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu May 21 15:39:23 EDT 2015

% figure of alternate GSL curves



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

ltrs='abcde';
ordr=[2 1 3 4 5];
wlabls={'ML_{2,1}','ML_{2,2}','ML_{1,1}','Gr','NC'};

clf; clear hp ht ht2;
for ii=1:length(ordr)
    iii=ordr(ii);
    hp(ii)=subplot(5,1,ii)
    wf=f2s{iii}(datsub,selmask);
    wsd=sd2s{iii}(datsub,selmask);

    
    clear plotdat;
    plotdat{1}.y = wf;
    plotdat{1}.dy = wsd;
    plotdat{1}.x = testX(datsub,3); 
    rgb=[0 0 0];
    
    [hl,hk]=PlotWithShadedErrors(plotdat,rgb,[],[],[],[-500 2010]);
    ylabel('GSL (mm)');
    box on;
    xlim([-500 2010]);
    ylim([-285 140]);
    yl=ylim;
    ht(ii)=text(-475,yl(end)-.1*diff(yl),[ltrs(ii) '    ' wlabls{ii}]);
    
end    

set(gcf,'PaperPosition',[.25 .25 8 10]);
set(ht,'fontweight','bold');
longticks(hp,2);

pdfwrite('GSL_diffpriors');