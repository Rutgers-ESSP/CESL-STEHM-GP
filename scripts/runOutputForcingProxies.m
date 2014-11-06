% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Nov 06 15:02:03 EST 2014


% plot of relevant forcing proxies


clear proxyT proxyO;

proxyT.time1 = Chesapeake.yr;
sub=find((Chesapeake.yr>=1000).*(Chesapeake.yr<=1800));
Chesapeake.offset=mean(Chesapeake.T(sub));
proxyT.Y = (Chesapeake.T-Chesapeake.offset)/10;
proxyT.dY = ones(size(Chesapeake.T))*.01;
proxyT.datid = ones(size(Chesapeake.T));
proxyT.sitelen = length(Chesapeake.T);

proxyT.time1 = [proxyT.time1 ; round(GOM.yr)];
sub=find((GOM.yr>=1000).*(GOM.yr<=1800));
GOM.offset=mean(GOM.T(sub));
proxyT.Y = [proxyT.Y ; (GOM.T-GOM.offset)/5];
proxyT.dY = [proxyT.dY ; ones(size(GOM.T))*.01];
proxyT.datid = [proxyT.datid ; 2*ones(size(GOM.T))];
proxyT.sitelen = [proxyT.sitelen ; length(GOM.T)];

proxyT.time2 = proxyT.time1; proxyT.meantime=proxyT.time1;
proxyT.limiting = zeros(size(proxyT.Y));
proxyT.compactcorr = zeros(size(proxyT.Y));
proxyT.lat = proxyT.datid; proxyT.long = proxyT.datid;
proxyT.Ycv = sparse(diag(proxyT.dY).^2);
proxyT.istg = zeros(size(proxyT.Y));
proxyT.siteid = [1 2]';
proxyT.sitenames={'Chesapeake','GOM'};
proxyT.sitecoords=[1 1 ; 2 2];

%%

sub=find(RothJoos.yr>-1000);
proxyO.time1 = [RothJoos.yr(sub)];
proxyO.Y = [(RothJoos.TSI(sub) - 1365)*1000];
proxyO.dY = [ ones(size(RothJoos.TSI(sub)))*1];
proxyO.datid = [ones(size(RothJoos.yr(sub)))*1];
proxyO.sitelen = [ length(sub)];

sub=find((haug.yr>-1000).*(~isnan(haug.Ti)));
proxyO.time1 = [proxyO.time1 ; round(haug.yr(sub))];
proxyO.Y = [proxyO.Y ; haug.Ti(sub)*1000];
proxyO.dY = [proxyO.dY ; ones(size(haug.Ti(sub)))*1];
proxyO.datid = [proxyO.datid ; ones(size(haug.Ti(sub)))*2];
proxyO.sitelen = [proxyO.sitelen ; length(haug.Ti(sub))];

proxyO.time2 = proxyO.time1; proxyO.meantime=proxyO.time1;
proxyO.limiting = zeros(size(proxyO.Y));
proxyO.compactcorr = zeros(size(proxyO.Y));
proxyO.lat = proxyO.datid; proxyO.long = proxyO.datid;
proxyO.Ycv = sparse(diag(proxyO.dY).^2);
proxyO.istg = zeros(size(proxyO.Y));
proxyO.siteid = [1 2]';
proxyO.sitenames={'TSI','Cariacao'};
proxyO.sitecoords=[1 1 ; 2 2];


smoothwindow=101;

[proxyTs,proxyTthet] = GPSmoothTideGauges(proxyT,smoothwindow,2.0,20)
[proxyOs,proxyOthet] = GPSmoothTideGauges(proxyO,smoothwindow,2.0,20)
[proxyTs] = GPSmoothTideGauges(proxyT,smoothwindow,0,20,proxyTthet(1,:));

sub=find((proxyTs.datid==1).*(proxyTs.meantime>=(min(Chesapeake.yr)+smoothwindow/2)).*(proxyTs.meantime<=(max(Chesapeake.yr)-smoothwindow/2)));
ChesapeakeS.yr=proxyTs.meantime(sub);
ChesapeakeS.Tanom=proxyTs.Y(sub)*10;
ChespeakeS.T=ChesapeakeS.Tanom+Chesapeake.offset;

sub=find((proxyTs.datid==2).*(proxyTs.meantime>=(min(GOM.yr)+smoothwindow/2)).*(proxyTs.meantime<=(max(GOM.yr)-smoothwindow/2)));
GOMS.yr=proxyTs.meantime(sub);
GOMS.Tanom=proxyTs.Y(sub)*5;
GOMS.T=GOMS.Tanom+GOM.offset;

sub=find((proxyOs.datid==1).*(proxyOs.meantime>=(-1000+smoothwindow/2)).*(proxyOs.meantime<=(max(RothJoos.yr)-smoothwindow/2)));
RothJoosS.yr=proxyOs.meantime(sub);
RothJoosS.TSI=proxyOs.Y(sub)/1000+1365;

sub=find((proxyOs.datid==2).*(proxyOs.meantime>=(-1000+smoothwindow/2)).*(proxyOs.meantime<=(max(haug.yr)-smoothwindow/2)));
haugS.yr=proxyOs.meantime(sub);
haugS.Ti=proxyOs.Y(sub)/1000;

colrs={'b','g','r'};
clear ha;
figure;
hp=subplot(2,1,1);
ha(1)=plot(ChesapeakeS.yr,ChesapeakeS.Tanom,colrs{1},'linew',2); hold on
ha(2)=plot(GOMS.yr,GOMS.Tanom,colrs{2},'linew',2); hold on
ha(3)=plot(-1000,-1000,colrs{3},'linew',2);
set(hp,'box','off');
ylim([min([ChesapeakeS.Tanom ; GOMS.Tanom])-.5  max([ChesapeakeS.Tanom ; GOMS.Tanom])+.5]);
xlabel('Year');
hyl=ylabel('Sea Surface Temperature (C)'); %set(hyl,'Color',colrs{1});
legend(ha,sprintf('Chesapeake - %0.1f C',Chesapeake.offset),sprintf('Gulf of Mexico - %0.1f C',GOM.offset),'Florida Current (Sv)','Location','Southwest')

hp(2)=axes('Position',get(hp(1),'Position'));
ha(3)=plot(Lundtransport.yr,Lundtransport.Sv,colrs{3},'linew',2); hold on
plot(Lundtransport.yr,Lundtransport.Sv+Lundtransport.dSv,[colrs{3} '--']); hold on
plot(Lundtransport.yr,Lundtransport.Sv-Lundtransport.dSv,[colrs{3} '--']); hold on
set(hp(2),'Color','none','XAxisLoc','top','YAxisLoc','right','box','off','xtickl',{''});
ylim([27.5 32]);
hyl=ylabel({'Transport (Sv)'}); set(hyl,'Color',colrs{3});

set(hp,'xlim',[-500 2000]);
pdfwrite('forcingproxies1');

figure;

hp=subplot(2,1,1);
plot(RothJoosS.yr,RothJoosS.TSI,colrs{1},'linew',2); hold on
set(hp,'box','off','ylim',[8 17]);
sub=find(RothJoosS.yr>=-500);
ylim([min(RothJoosS.TSI(sub))-.05 max(RothJoosS.TSI(sub))+.05]);
hyl=ylabel('TSI (W/m^2)'); set(hyl,'Color',colrs{1});
xlabel('Year');

hp(2)=axes('Position',get(hp(1),'Position'));
plot(haugS.yr,haugS.Ti,colrs{2},'linew',1); hold on
set(hp(2),'Color','none','XAxisLoc','top','YAxisLoc','right','box','off','xtickl',{''});
sub=find(haugS.yr>=-500);
ylim([min(haugS.Ti(sub))-.02  max(haugS.Ti(sub))+.02]);
hyl=ylabel({'Cariaco Basin % Ti','(lower = dryer)'}); set(hyl,'Color',colrs{2});

set(hp,'xlim',[-500 2000]);
pdfwrite('forcingproxies2');

M=[ChesapeakeS.yr,ChesapeakeS.T];
save -ascii ChesapeakeS.asc M

M=[GOMS.yr,GOMS.T];
save -ascii GOM.asc M

M=[RothJoosS.yr,RothJoosS.TSI];
save -ascii RothJoosS.asc M

M=[haugS.yr,haugS.Ti];
save -ascii haugS.asc M
