function [hp,outtable]=PlotPaleoRSLOneSite(testt,fs,sds,Vs,difftimestep,PX,YdD,sitename,thetL)

% [hp,outtable]=PlotPaleoRSLOneSite(testt,fs,sds,[Vs],[difftimestep],PX,YdD,[sitename],[thetL])
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Apr 13 18:38:07 EDT 2014

defval('scaler',1e-3);
defval('do2s',1);
defval('difftimestep',0);
defval('sitename','');
defval('thetL',[]);
defval('Vs',diag(sds.^2));

datid=PX.datid;
istg=PX.istg;
limiting=PX.limiting;
compactcorr=PX.compactcorr;
time1=PX.time1;
time2=PX.time2;
dY=PX.dY;
meantime=PX.meantime;
defval('YdD',PX.Y);


if difftimestep>0
	Mdiff = bsxfun(@eq,testt,testt')-bsxfun(@eq,testt,testt'+difftimestep);
	sub=find(sum(Mdiff,2)==0);
	Mdiff=Mdiff(sub,:);
	difftimes=bsxfun(@rdivide,abs(Mdiff)*testt,sum(abs(Mdiff),2));;
	Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testt);

	dfs=Mdiff*fs;
	dVs=Mdiff*Vs*Mdiff';
	dsds=sqrt(diag(dVs));
end

clf;
hp(1)=subplot(2,1,1)
hold on;

plot(testt,scaler*(fs),'r-','linew',2); hold on;
plot(testt,scaler*(fs+sds),'r--'); 
plot(testt,scaler*(fs-sds),'r--'); 
if do2s
	plot(testt,scaler*(fs+2*sds),'r:'); 
	plot(testt,scaler*(fs-2*sds),'r:');
end


trainsub0=find((limiting==0).*(istg<=1).*(compactcorr~=0));
for jj=trainsub0'
%	plot([time1(jj) time2(jj)],scaler*[1 1]*YdD(jj),'g'); hold on;
%	plot(meantime(jj)*[1 1],scaler*(YdD(jj)+2*[-1 1]*dY(jj)),'g');
	plot(meantime(jj)*[1 1],scaler*(YdD(jj)+2*[-1 1]*sqrt(dY(jj).^2+(YdD(jj)*thetL(end)).^2)),'g');
end
%plot(meantime(trainsub0),scaler*YdD(trainsub0),'g.');

%trainsub0=find((limiting==0).*(istg<=1).*(compactcorr==0));
trainsub0=find((limiting==0).*(istg<=1));
for jj=trainsub0'
	plot([time1(jj) time2(jj)],scaler*[1 1]*YdD(jj),'b'); hold on;
	plot(meantime(jj)*[1 1],scaler*(YdD(jj)+2*[-1 1]*dY(jj)),'b');
end
plot(meantime(trainsub0),scaler*YdD(trainsub0),'b.');



trainsub0=find((limiting==1));
for jj=trainsub0'
	plot([time1(jj) time2(jj)],scaler*[1 1]*YdD(jj),'c');
	plot(meantime(jj)*[1 1],scaler*(YdD(jj)+2*[-1 1]*dY(jj)),'c');
end
plot(meantime(trainsub0),scaler*YdD(trainsub0),'cv');

trainsub0=find((limiting==-1));
for jj=trainsub0'
	plot([time1(jj) time2(jj)],scaler*[1 1]*YdD(jj),'m');
	plot(meantime(jj)*[1 1],scaler*(YdD(jj)+2*[-1 1]*dY(jj)),'m');
end
plot(meantime(trainsub0),scaler*YdD(trainsub0),'m^');

xlabel('Year CE'); ylabel('m');
xlim([min(time1) max(testt)]);
ylim([min(union(floor(scaler*(YdD-2*dY)),floor(scaler*(fs-2*sds)))) max(0,max(union(ceil(scaler*(YdD+2*dY)),ceil(scaler*(fs+2*sds))))) ]);
box on;


hp(2)=subplot(2,1,2);
if difftimestep>0
	scaler = 1;

	plot([-10000 10000],[0 0],'k-'); hold on;
	plot(difftimes,scaler*(dfs),'r-','linew',2); hold on;
	plot(difftimes,scaler*(dfs+dsds),'r--'); 
	plot(difftimes,scaler*(dfs-dsds),'r--'); 
	plot(difftimes,scaler*(dfs+2*dsds),'r:'); 
	plot(difftimes,scaler*(dfs-2*dsds),'r:');
	xlabel('Year CE'); ylabel(['mm/y (' num2str(difftimestep) '-y avg.)']);

	ylim([min(0,min(floor(scaler*dfs-dsds))) max(ceil(scaler*dfs+dsds))]);
end

set(hp(1:2),'xlim',[min(time1) max(testt)]);
hp(3)=axes('position',get(hp(1),'position'));
hp(4)=axes('position',get(hp(2),'position'));
set(hp(3:4),'color','none','xaxisloc','top','yaxisloc','right');
set(hp(1:2),'box','off');
set(hp(3:4),'ylim',get(hp(1),'ylim')); set(hp(3:4),'xdir','reverse','xlim',1950-[max(testt) min(time1)]);
xticks=get(hp(1),'xtick');
set(hp(3:4),'xtick',1950-xticks(end:-1:1));

set(hp(3:4),'yticklabel',{}); xlabel(hp(3),'Year BP'); xlabel(hp(4),'Year BP');
longticks(hp,2);

if difftimestep>0
	set(hp([1 4]),'xticklabel',{}); xlabel(hp(1),''); xlabel(hp(4),'');
	movev(hp([2 4]),.09);
else
	delete(hp([2 4]));
end

outtable=[sitename '\n'];
outtable=[outtable ['Gaussian process analysis\n' datestr(now) '\n\n']];
if length(thetL)>0
	outtable=[outtable 'theta:,'];
	outtable=[outtable sprintf('%3.5f,',thetL)];
end
outtable=[outtable '\n\nYears CE,Years BP,Sea level (mm),+/- 2 sigma\n'];
u1=[testt(:) 1950-testt(:) fs(:) 2*sds(:)]';
outtable=[outtable sprintf('%6.0f,%6.3f,%6.3f,%6.3f\n',u1)];
outtable=[outtable '\n\n'];

if difftimestep>0
	outtable=[outtable '\n\nYears CE,Sea level rise (mm/y) [' num2str(difftimestep) '-y avg],+/- 2 sigma\n'];
	u1=[difftimes(:) dfs(:) 2*dsds(:)]';
	outtable=[outtable sprintf('%6.0f,%6.3f,%6.3f\n',u1)];
end
