% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 18 00:04:18 EST 2014


addpath('~/Dropbox/Code/TGAnalysis/MFILES');
addpath([pwd '/MFILES']);
IFILES=[pwd '/IFILES'];
addpath(pwd)

WORKDIR='140217';
if ~exist(WORKDIR,'dir')
	mkdir(WORKDIR);
end
cd(WORKDIR);


%
runImportHoloceneDataSets;
runSetupHoloceneCovariance;

trainsets=[1:3]; % train w/TG+GSL+PX, w/TG+PX, w/PX
trainspecs=ones(size(trainsets));

clear thetTGG trainsubsubset;
for ii=1:length(trainsets)
    [thetTGG{ii},trainsubsubset{ii}]=OptimizeHoloceneCovariance(datasets{trainsets(ii)},modelspec(trainspecs(ii)));
end

%%%%%%%%

fid=fopen('thetTGG.tsv','w');
for ii=1:length(datasets)
    fprintf(fid,[datasets{ii}.label '\t']);
    fprintf(fid,'%0.3f\t',thetTGG{ii});
    fprintf(fid,'\n');
end
fclose(fid);


%%%%

fid=fopen('TGandProxyData.tsv','w');
fprintf(fid,'ID\t1ime1\ttime2\tlimiting\tY-GIAproj\tY\tdY\tcompactcorr\tistg\tlat\tlong\n');
ii=1;
for i=1:size(datasets{ii}.datid,1)
    compactcorr=full(datasets{ii}.compactcorr);
	fprintf(fid,'%d\t',datasets{ii}.datid(i));
	fprintf(fid,'%0.1f\t',datasets{ii}.time1(i));
	fprintf(fid,'%0.1f\t',datasets{ii}.time2(i));
	fprintf(fid,'%d\t',datasets{ii}.limiting(i));
	fprintf(fid,'%0.2f\t',datasets{ii}.Y(i));
	fprintf(fid,'%0.2f\t',datasets{ii}.Y(i)+datasets{ii}.GIAproj(i));
	fprintf(fid,'%0.2f\t',datasets{ii}.dY(i));
	fprintf(fid,'%0.2f\t',compactcorr(i));
	fprintf(fid,'%d\t',datasets{ii}.istg(i));
	fprintf(fid,'%0.2f\t',datasets{ii}.lat(i));
	fprintf(fid,'%0.2f\n',datasets{ii}.long(i));
end
fclose(fid)

for ii=1:length(datasets)
    datasets{ii}.dY = sqrt(datasets{ii}.dY0.^2 + (thetTGG{ii}(end)*datasets{ii}.compactcorr).^2);
    datasets{ii}.Ycv = datasets{ii}.Ycv0 + diag(thetTGG{ii}(end)*datasets{ii}.compactcorr).^2;
end

%% now do a regression

testsitedef.sites=[
0   1e6   1e6
idFL latFL longFL
idNC+1 latNC1 longNC1
idNC latNC longNC
idNC+2 latNC2 longNC2
idNJ+1 latNJ1 longNJ1
idNJ latNJ longNJ
idNJ+2 latNJ2 longNJ2
idNS latNS longNS
];
testsitedef.names={'GSL','FL','NC_TP','NC','NC_SP','NJ_CM','NJ','NJ_LP','NS'};
testsitedef.names2={'GSL','Nassau','Tump Point','North Carolina','Sand Point','Cape May','New Jersey','Leeds Point','Chezzetcook'};

noiseMasks = ones(8,length(thetTGG{1}));
noiseMasks(2,[8 10] )=0; %without linear
noiseMasks(3,[1 12])=0; %only regional and local
noiseMasks(4,[1 12 8 10])=0; %only regional and local non-linear
noiseMasks(5,[1 2 3 4])=0; %only regional and local linear
noiseMasks(6,[1 3 4 8 10 13])=0; %Greenland only
noiseMasks(7,[1 12 11 4 10])=0; %only regional
noiseMasks(8,[1 12 11 4 10 8])=0; %only regional non-linear
noiseMasklabels={'full','nonlin','regloc','reglocnonlin','regloclin','gis','reg','regnonlin'};

GISfpt.lat=GISfplat;
GISfpt.long=GISfplong;
GISfpt.fp=GISfp;

ICE5G.lat=ICE5Glat;
ICE5G.long=ICE5Glon;
ICE5G.gia=ICE5Ggia;


testt = [-2000:10:2010];

for ii=1:length(datasets)
for jj=1:length(trainsets)

    labls{ii,jj}=['_' datasets{ii}.label '_tr' datasets{trainsets(jj)}.label];
    disp(labls{ii,jj});

    trainsub = find((datasets{ii}.limiting==0)); % only index points
    [f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj}]=RegressHoloceneDataSets(datasets{ii},testsitedef,modelspec(trainspecs(jj)),thetTGG{jj},GISfpt,ICE5G,trainsub,trainsubsubset{ii},noiseMasks,testt,refyear)

end
end

save ~/tmp/CESL

%%

for ii=1:length(datasets)
for jj=1:length(trainsets)

    labl=labls{ii,jj}; disp(labl);
   makeplots_sldecomp(datasets{ii},f2s{ii,jj},sd2s{ii,jj},V2s{ii,jj},testlocs{ii,jj},labl);
    
    testreg=testlocs{ii,jj}.reg;
    testsites=testlocs{ii,jj}.sites;
    testX=testlocs{ii,jj}.X;
    testnames2=testlocs{ii,jj}.names2;

    firstyears=[800 0 -500 1800 1900];
    lastyears=[1800 1800 1700 1900 2000];

    for kk=1:size(testsites,1)
        for pp=1:length(firstyears)
            sub=find((testreg==testsites(kk,1)).*(testX(:,3)>=firstyears(pp)).*(testX(:,3)<=lastyears(pp)));
            M = zeros(1,length(sub)); M(1)=-1; M(end)=1; 
            M=M/(testX(sub(end),3)-testX(sub(1),3));

            if pp==1
                fslopelin(kk,ii,jj) = M*f2s{ii,jj}(sub,5);
                Vslope = M*V2s{ii,jj}(sub,sub,5)*M';
                sdslopelin(kk,ii,jj)=sqrt(diag(Vslope));
            end

            fslopeavg(kk,ii,jj,pp) = M*f2s{ii,jj}(sub,1);
            Vslope = M*V2s{ii,jj}(sub,sub,1)*M';
            sdslopeavg(kk,ii,jj,pp)=sqrt(diag(Vslope));

            fslopeavglessGSL(kk,ii,jj,pp) = M*f2s{ii,jj}(sub,3);
            Vslope = M*V2s{ii,jj}(sub,sub,3)*M';
            sdslopeavglessGSL(kk,ii,jj,pp)=sqrt(diag(Vslope));

       end

    end
    
    
    fid=fopen(['linrates' labl '.tsv'],'w');
    fprintf(fid,['Regional+Local Linear Rates (mm/y), ' labl '\n']);
    fprintf(fid,'Site\tRate (linear)\t2s');
    for pp=1:length(firstyears)
        fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
    end
    for pp=1:length(firstyears)
        fprintf(fid,'\tRate (avg w/o GSL, %0.0f-%0.0f)\t2s',[firstyears(pp) lastyears(pp)]);
    end
    fprintf(fid,'\n');
    for kk=1:size(testsites,1)
        fprintf(fid,testnames2{kk});
        fprintf(fid,'\t%0.2f',[fslopelin(kk,ii,jj) 2*sdslopelin(kk,ii,jj)]);
        for pp=1:length(firstyears)
            fprintf(fid,'\t%0.2f',[fslopeavg(kk,ii,jj,pp) 2*sdslopeavg(kk,ii,jj,pp)]);
        end
        for pp=1:length(firstyears)
            fprintf(fid,'\t%0.2f',[fslopeavglessGSL(kk,ii,jj,pp) 2*sdslopeavglessGSL(kk,ii,jj,pp)]);
        end
        fprintf(fid,'\n');
    end
    fclose(fid);

        
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% slopes plot -- TBD

fslope=fslope(:);
sub=find(testsites(:,2)<1e3);
sub=setdiff(sub,[4 7]);

figure;
% plot observation sites
scatter((testsites(sub,3)),testsites(sub,2),100,fslope(sub),'filled'); hold on;

xl=xlim;
yl=ylim;

clf;
usstates = shaperead('usastatehi', 'UseGeoCoords', true);
plot(([usstates.Lon]),[usstates.Lat],'k','linew',.5); hold on;

% Canadian province boundaries from http://www.nws.noaa.gov/geodata/catalog/national/html/province.htm

if exist([IFILES '/province/PROVINCE.SHP']); canada = shaperead([IFILES '/province/PROVINCE.SHP'],'UseGeoCoords',true); else; canada=[]; end
if length(canada)>0 ; plot(([canada.Lon]),[canada.Lat],'k','linew',.5); end

% plot tide gauge locations
trainsub2=find((istg==1).*(lat<1000));
[u,ui]=unique(datid(trainsub2));
scatter(long(trainsub2(ui)),lat(trainsub2(ui)),30,'k','filled'); hold on;

% plot observation sites
scatter((testsites(sub,3)),testsites(sub,2),100,fslope(sub),'filled'); hold on;


%[~,handl]=plotcont([xl(1) yl(2)],[xl(2) yl(1)],1); hold on; set(handl,'linew',.5);
axis tight;
xlim(xl);
ylim(yl);
colorbar;
box on;
%title('Sea level anomaly rates, 900-1900 CE (mm/y)');
pdfwrite('map_linrates');



%%% 

sitesub=[2 3 5 6 8 9];
starttimes=[];
endtimes=[];
colrs={'r','c','b','g','m','y'};

starttimes=[1:size(testsites,1)]*-Inf;
endtimes=[1:size(testsites,1)]*Inf;

for i=sitesub
	subA = find(testreg == testsites(i,1));
	distfrom=dDist(testsites(i,2:3),[lat long]);
	subB=find(distfrom<.1);
	starttimes(i) = max(min(testt),min(time1(subB)));
	endtimes(i) = max(time2(subB));
end

starttimes(4)=min(starttimes([3 5]));
starttimes(7)=min(starttimes([6 8]));


datsub=find(ismember(testreg,testsites(sitesub,1)));

% figure of sites with data, full curves
% figure of sites with data, without linear
% figure of sites with data, without linear or GSL

for selmask=[1 2 3 4]
	
	clf;
	
	wf=Mref*f2s(:,selmask);
	wV=Mref*V2s(:,:,selmask)*Mref';
	wsd=sqrt(diag(wV));
	[hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wsd(datsub),colrs,starttimes(sitesub),endtimes(sitesub),0,160);
	
	if selmask==1
		legend(hl,testnames2{sitesub},'Location','Southeast');
	else
		legend(hl,testnames2{sitesub},'Location','Southwest');
	end
	set(hp,'xlim',[-1000 2010]);
%	set(hp(2),'ylim',[-2 5]);

	pdfwrite(['paleorsl_' noiseMasklabels{selmask}]);

	fid=fopen(['paleorsl_' noiseMasklabels{selmask} '.tsv'],'w');
	fprintf(fid,outtable);
	fclose(fid);
end


% figure of regional averages, without linear or GSL

sitesub = [2 4 7 9];
colrs={'r','b','g','y'};
selmask=8;


datsub=find(ismember(testreg,testsites(sitesub,1)));

%wf=f2s(:,selmask);
%wsd=sd2s(:,selmask);

wf=Mref*f2s(:,selmask);
wV=Mref*V2s(:,:,selmask)*Mref';
wsd = sqrt(diag(wV));

clf;
[hp,hl,~,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wsd(datsub),colrs,starttimes(sitesub),endtimes(sitesub),0,40);
legend(hl,testnames2{sitesub},'Location','Southwest');

set(hp,'xlim',[-1000 2010]);

pdfwrite(['paleorsl_' noiseMasklabels{selmask}]);

fid=fopen(['paleorsl_' noiseMasklabels{selmask} '.tsv'],'w');
fprintf(fid,outtable);
fclose(fid);

% figure of GSL and rate of GSL change

sitesub=[1];
colrs={'k'};
selmask=1;

wf=Mref*f2s(:,selmask);
wV=Mref*V2s(:,:,selmask)*Mref';
wsd=sqrt(diag(wV));

datsub=find(ismember(testreg,testsites(sitesub,1)));

for timesteps=[100 1000 160 40 20 2200]

    clf;
    [hp,hl,hl2,dGSL,dGSLsd,dGSLV,outtable,difftimes,diffreg]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,starttimes(sitesub),endtimes(sitesub),0,timesteps);
    set(hp,'xlim',[-1000 2010]);

    pdfwrite(['GSL_' num2str(timesteps) 'y']);


    fid=fopen(['GSL_' num2str(timesteps) 'y.tsv'],'w');
    fprintf(fid,outtable);

    Nsamps=10000;
    sub1=find((difftimes<=1800).*(difftimes>=1500)); sub2=find(difftimes>1800);
    sub1a=find((difftimes<=1800).*(difftimes>=0)); sub2=find(difftimes>1800);   
    samps=mvnrnd(dGSL,dGSLV,Nsamps);
    if (length(sub1)>0).*(length(sub2)>0)
        [max1,max1i]=max(samps(:,sub1),[],2);
        [max1a,max1ia]=max(samps(:,sub1a),[],2);
        difr=bsxfun(@minus,samps(:,sub2),max1);
        difra=bsxfun(@minus,samps(:,sub2),max1a);
        Ppositive = sum(samps(:,sub2)>0,1)/size(samps,1);
        q=cumprod((samps(:,sub2(end:-1:1))>0)*1,2); q=q(:,end:-1:1);
        Ppositiveseries = sum(q,1)/size(samps,1);
        Plarger = sum(difr>0,1)/size(samps,1);
        Plargera = sum(difra>0,1)/size(samps,1);
        fprintf(fid,'\n\nProbability faster than during fastest interval\n');
        fprintf(fid,'Central year\tP > 0\tP all subsequent > 0\tP centered 1500-1800\tP centered 0-1800');
        fprintf(fid,'\n%0.0f\t%0.3f\t%0.3f\t%0.3f\t%0.3f',[difftimes(sub2)' ; Ppositive; Ppositiveseries ; Plarger ; Plargera]);
    end

    fclose(fid);
    
end

% sample GSL rate curve

Nsamps=10000;
suboverallyrs=find(mod(difftimes,100)==50);
samps=mvnrnd(dGSL(suboverallyrs),dGSLV(suboverallyrs,suboverallyrs),Nsamps);
subyrs=find(difftimes(suboverallyrs)>=1800);
clear lastyrlarger;
for ii=1:length(subyrs)
	curyr=difftimes(suboverallyrs(subyrs(ii)));
	for jj=1:Nsamps
		sub=find(difftimes(suboverallyrs)<curyr);
		sub2=find(samps(jj,sub)>samps(jj,sub(end)+1));
		if length(sub2)>0
			lastyrlarger(jj,ii)=difftimes(suboverallyrs(sub2(end)));
		else
			lastyrlarger(jj,ii)=min(difftimes)-1;
		end		
	end
end
quantile(lastyrlarger,.95)

% sample GSL rate curve

Nsamps=10000;
samps=mvnrnd(wf(datsub),wV(datsub,datsub),Nsamps);
subyrs=find(testX(datsub,3)>=1800);
clear lastyrlarger;
for ii=1:length(subyrs)
	curyr=testX(datsub(subyrs(ii)),3);
	for jj=1:Nsamps
		sub=find(testX(datsub,3)<curyr);
		sub2=find(samps(jj,sub)>samps(jj,sub(end)+1));
		if length(sub2)>0
			lastyrlarger(jj,ii)=testX(datsub(sub2(end)),3);
		else
			lastyrlarger(jj,ii)=min(testX(datsub,3))-1;
		end		
	end
end
quantile(lastyrlarger,.95)


%%%

% figure of GIS and rate of GIS change

sitesub=[1];
colrs={'k'};
selmask=6;

wf=Mref*f2s(:,selmask);
wV=Mref*V2s(:,:,selmask)*Mref';
wsd=sqrt(diag(wV));

datsub=find(ismember(testreg,testsites(sitesub,1)));

clf;
[hp,hl,hl2,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,[],[],0,100);
set(hp,'xlim',[-1000 2010]);

pdfwrite(['GIS']);

fid=fopen('GIS.tsv','w');
fprintf(fid,outtable);
fclose(fid);

clf;
[hp,hl,hl2,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),wf(datsub),wV(datsub,datsub),colrs,[],[],0,600);
set(hp,'xlim',[-1000 2010]);

pdfwrite(['GIS_600y']);

fid=fopen('GIS_600y.tsv','w');
fprintf(fid,outtable);
fclose(fid);

%%%%

% figure of sea-level differences


clear subs;
for i=1:size(testsites,1)
	subs{i}=find(testreg==testsites(i,1));
end
Mgrad=zeros(length(subs{1}),length(testX));
for i=1:size(testsites,1)
	Mgrad(subs{i},subs{i}) = -eye(length(subs{i}));
	Mgrad(subs{i},subs{7}) = Mgrad(subs{i},subs{7}) + eye(length(subs{7}));
	gradstarttimes(i) = max(starttimes([i 7]));
end

gradf = Mgrad*f2s(:,2);
gradV = Mgrad*V2s(:,:,2)*Mgrad';

gradf=Mref*gradf;
gradV=Mref*gradV*Mref';
gradsd=sqrt(diag(gradV));

sitesub=[2 4 9];
datsub=find(ismember(testreg,testsites(sitesub,1)));
colrs={'r','g','y'};

clf;
[hp,hl,hl2,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),gradf(datsub),gradV(datsub,datsub),colrs,gradstarttimes(sitesub),[],0,40);
set(hp,'xlim',[-1000 2010]);
legend(hl,'FL-NJ','NC-NJ','NS-NJ','Location','Southwest');

pdfwrite(['RSLgrad']);

fid=fopen('RSLgrad.tsv','w');
fprintf(fid,outtable);
fclose(fid);

clf;
[hp,hl,hl2,~,~,~,outtable]=PlotPSLOverlay(testX(datsub,3),testreg(datsub),testsites(sitesub,1),gradf(datsub),gradV(datsub,datsub),colrs,gradstarttimes(sitesub),[],0,600);
set(hp,'xlim',[-1000 2010]);
pdfwrite(['RSLgrad_600y']);

fid=fopen('RSLgrad_600y.tsv','w');
fprintf(fid,outtable);
fclose(fid);

%%%


%%



MlessGSL=eye(length(testreg),length(testreg));
for i=1:size(testsites,1)
	MlessGSL(subs{i},find(testreg==0)) = MlessGSL(subs{i},find(testreg==0))-eye(length(subs{1}));
end

clear gradfs gradVs gradsds;
for i=1:size(noiseMasks,1)
	gradfs(:,i) = Mgrad*f2s(:,i);
	gradVs(:,:,i) = Mgrad * V2s(:,:,i) * Mgrad';
	gradsds(:,i) = sqrt(diag(gradVs(:,:,i)));
end

gradtime = abs(Mgrad)*testX(:,3)/2;

clear flessGSL VlessGSL sdlessGSL;
for i=1:size(noiseMasks,1)
	flessGSL(:,i) = MlessGSL*f2s(:,i);
	VlessGSL(:,:,i) = MlessGSL * V2s(:,:,i) * MlessGSL';
	sdlessGSL(:,i) = sqrt(diag(VlessGSL(:,:,i)));
end


%

% GSL minimum

Nsamps=1000;
sub=find(testreg==0);
samps=mvnrnd(f2s(sub,1),V2s(sub,sub,1),Nsamps);
sub2=find(testt>1600);
[m,mi]=min(samps(:,sub2),[],2);



for difftimestep=[160 40]
	timestr = ['_' num2str(difftimestep) 'y'];

	Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+difftimestep);
	Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
	sub=find(sum(Mdiff,2)==0);
	Mdiff=Mdiff(sub,:);
	difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));;
	diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));;
	Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));
	df1=Mdiff*f1;
	dV1=Mdiff*V1*Mdiff';
	dsd1=sqrt(diag(dV1));

	clear df2s dV2s dsd2s;
	for i=1:size(noiseMasks,1)
		df2s(:,i)=Mdiff*f2s(:,i);
		dV2s(:,:,i)=Mdiff*V2s(:,:,i)*Mdiff';
		dsd2s(:,i)=sqrt(diag(dV2s(:,:,i)));
	end


	Mdiff2=Mdiff(2:length(gradfs)-1,1:length(gradfs));
	dgradtime=bsxfun(@rdivide,abs(Mdiff2)*gradtime,sum(abs(Mdiff2),2));;
	sub=find(~isnan(dgradtime));

	Mdiff2 = Mdiff2(sub,:);
	dgradtime=dgradtime(sub,:);
	clear dgradfs dgradVs dgradsds;
	for i=1:size(noiseMasks,1)
		dgradfs(:,i)=Mdiff2*gradfs(:,i);
		dgradVs(:,:,i)=Mdiff2*gradVs(:,:,i)*Mdiff2';
		dgradsds(:,i)=sqrt(diag(dgradVs(:,:,i)));
	end

	for j=1:size(testsites,1)
		clear hp;
		dists = sqrt((lat-testsites(j,2)).^2 + (long-testsites(j,3)).^2);
		trainsub2 = find((dists<1).*(limiting==0));
		trainsub2=intersect(trainsub2,trainsubsubset);
		testsub = find(testreg==testsites(j,1));
		difftestsub = find(diffreg==testsites(j,1));
		figure;
		hp(1)=subplot(2,1,1)
		for i=trainsub2(:)'
			plot([time1(i) time2(i)],[1 1]*Y(i)/1000);
			hold on;
			plot([1 1]*meantime(i),(Y(i)+2*dY0(i)*[-1 1])/1000);
		end
		if length(trainsub2)>0
			mintime = min(time1(trainsub2));
		end

		hpl(2)=	plot(testX(testsub,3),1e-3*(f2s(testsub,1)),'g'); hold on;
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)+sd2s(testsub,1)),'g--');
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)-sd2s(testsub,1)),'g--');
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)+2*sd2s(testsub,1)),'g:');
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)-2*sd2s(testsub,1)),'g:');
		xlim([-1000 2020]);
		ylabel('m');
		title(testnames{j});
	%	legend(hpl,'w/o age unc','w/ age unc','Location','northwest');

		hp(2)=subplot(2,1,2)

		plot(difftimes(difftestsub),df2s(difftestsub,1),'g','linew',2); hold on
		plot(difftimes(difftestsub),df2s(difftestsub,1)+dsd2s(difftestsub,1),'g--');
		plot(difftimes(difftestsub),df2s(difftestsub,1)-dsd2s(difftestsub,1),'g--');
		plot(difftimes(difftestsub),df2s(difftestsub,1)+2*dsd2s(difftestsub,1),'g:');
		plot(difftimes(difftestsub),df2s(difftestsub,1)-2*dsd2s(difftestsub,1),'g:');
		ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
		xlabel('Year CE');

		set(hp,'xlim',[-1000 2020]);
		pdfwrite([testnames{j} timestr]);

	end

	subsiteslist={[5 2 1 3]};
	subsitesname={'NSNJNCFL'};

	for nn=1:length(subsiteslist)
		clear hpl;
		subsites=subsiteslist{nn};
		colrs='rgbkcmy';
		figure;
		hp(2)=subplot(3,1,2);
		for jj=1:length(subsites)
			i=subsites(jj);
			if testsites(i,1)>0
				difftestsub = find(diffreg==testsites(i,1));
				difftestsub = intersect(difftestsub,find(difftimes>=minagesite(i)));

				hpl(jj)=plot(difftimes(difftestsub),df2s(difftestsub,1),[colrs(jj)],'linew',2); hold on
				plot(difftimes(difftestsub),df2s(difftestsub,1)+dsd2s(difftestsub,1),[colrs(jj) '--']);
				plot(difftimes(difftestsub),df2s(difftestsub,1)-dsd2s(difftestsub,1),[colrs(jj) '--']);
				plot(difftimes(difftestsub),df2s(difftestsub,1)+2*dsd2s(difftestsub,1),[colrs(jj) ':']);
				plot(difftimes(difftestsub),df2s(difftestsub,1)-2*dsd2s(difftestsub,1),[colrs(jj) ':']);
			end
		end
		title('Local sea level');
		ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
		xlabel('Year CE');

		hp(1)=subplot(3,1,1);

		jj=jj+1;
		difftestsub = find(diffreg==0);
		hpl(jj)=plot(difftimes(difftestsub),df2s(difftestsub,1),[colrs(jj)],'linew',2); hold on
		plot(difftimes(difftestsub),df2s(difftestsub,1)+dsd2s(difftestsub,1),[colrs(jj) '--']);
		plot(difftimes(difftestsub),df2s(difftestsub,1)-dsd2s(difftestsub,1),[colrs(jj) '--']);
		plot(difftimes(difftestsub),df2s(difftestsub,1)+2*dsd2s(difftestsub,1),[colrs(jj) ':']);
		plot(difftimes(difftestsub),df2s(difftestsub,1)-2*dsd2s(difftestsub,1),[colrs(jj) ':']);

		ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
		xlabel('Year CE');
		legend(hpl,testnames{[subsites 4]},'Location','Northwest');
		title('Global sea level');

		hp(3)=subplot(3,1,3);
		for jj=1:length(subsites)
			i=subsites(jj);
			if testsites(i,1)>0
				difftestsub = find(diffreg==testsites(i,1));
				difftestsub = intersect(difftestsub,find(difftimes>=minagesite(i)));
				
				hpl(jj)=plot(difftimes(difftestsub),df2s(difftestsub,3),[colrs(jj)],'linew',2); hold on
				plot(difftimes(difftestsub),df2s(difftestsub,3)+dsd2s(difftestsub,3),[colrs(jj) '--']);
				plot(difftimes(difftestsub),df2s(difftestsub,3)-dsd2s(difftestsub,3),[colrs(jj) '--']);
				plot(difftimes(difftestsub),df2s(difftestsub,3)+2*dsd2s(difftestsub,3),[colrs(jj) ':']);
				plot(difftimes(difftestsub),df2s(difftestsub,3)-2*dsd2s(difftestsub,3),[colrs(jj) ':']);
			end
		end

		title('Local sea level anomaly');
		ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
		xlabel('Year CE');

		set(hp,'xlim',[-1000 2020]);
		pdfwrite([subsitesname{nn} 'rate' timestr]);


	end
	

	%

	sub=find(floor(datid/1000)==floor(idNJ/1000));
	minageNJNC=min(time1(sub));
	sub=find(floor(datid/1000)==floor(idNC/1000));
	minageNJNC=max(minageNJNC,min(time1(sub)));
	
	figure; clear hp;
	hp(1)=subplot(2,1,1);
	sub=find(testt>=minageNJNC);
	plot(testt(sub),(gradfs(sub,1))/1000,'g','linew',2); hold on;
	plot(testt(sub),(gradfs(sub,1)+gradsds(sub,1))/1000,'g--');
	plot(testt(sub),(gradfs(sub,1)-gradsds(sub,1))/1000,'g--');
	plot(testt(sub),(gradfs(sub,1)+2*gradsds(sub,1))/1000,'g:');
	plot(testt(sub),(gradfs(sub,1)-2*gradsds(sub,1))/1000,'g:');
	title('NJ - NC'); ylabel('m');

	hp(2)=subplot(2,1,2);
	sub=find(dgradtime>=minageNJNC);
	plot(dgradtime(sub),dgradfs(sub,1),'g','linew',2); hold on;
	plot(dgradtime(sub),dgradfs(sub,1)+dgradsds(sub,1),'g--');
	plot(dgradtime(sub),dgradfs(sub,1)-dgradsds(sub,1),'g--');
	plot(dgradtime(sub),dgradfs(sub,1)+2*dgradsds(sub,1),'g:');
	plot(dgradtime(sub),dgradfs(sub,1)-2*dgradsds(sub,1),'g:');
	title('NJ - NC');
	ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);


	set(hp,'xlim',[0 2020]);
	pdfwrite(['NCNJgrad' timestr]);



	clf; clear hp hpl;
	hp(1)=subplot(3,1,1);

	i=1;
	difftestsub = find(diffreg==0);
	hpl(i)=plot(difftimes(difftestsub),df2s(difftestsub,1),[colrs(i)],'linew',2); hold on
	plot(difftimes(difftestsub),df2s(difftestsub,1)+dsd2s(difftestsub,1),[colrs(i) '--']);
	plot(difftimes(difftestsub),df2s(difftestsub,1)-dsd2s(difftestsub,1),[colrs(i) '--']);
	plot(difftimes(difftestsub),df2s(difftestsub,1)+2*dsd2s(difftestsub,1),[colrs(i) ':']);
	plot(difftimes(difftestsub),df2s(difftestsub,1)-2*dsd2s(difftestsub,1),[colrs(i) ':']);



	ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
	title('Global sea level');


	Mavg=sparse(abs(bsxfun(@minus,Marcottgl_yr,Marcottgl_yr'))<difftimestep/2);
	sub=find(sum(Mavg,2)==max(sum(Mavg,2)));
	Mavg=Mavg(sub,:)/max(sum(Mavg,2));
	Marcottgl_avgT = Mavg*Marcottgl_T;
	Marcottgl_avgyr = Mavg*Marcottgl_yr;
	%Marcottgl_avgdT = sqrt(diag(Mavg*diag(Marcottgl_dT.^2)*Mavg'));
	Marcottgl_avgdT = Marcottgl_dT(sub);

	Mavg=sparse(abs(bsxfun(@minus,Mann_yr,Mann_yr'))<difftimestep/2);
	sub=find(sum(Mavg,2)==max(sum(Mavg,2)));
	Mavg=Mavg(sub,:)/max(sum(Mavg,2));
	Mann_avgT = Mavg*Mann_T;
	Mann_avgyr = Mavg*Mann_yr;




	hp(2)=axes('position',get(hp(1),'position'));

	i=2;
	hpl(i)=plot(Mann_avgyr,Mann_avgT,[colrs(i)]); hold on;


	%i=3;
	%hpl(i)=plot(Marcottgl_yr,Marcott_Next_T,[colrs(i)]); hold on;
	%plot(Marcottgl_yr,Marcott_Next_T+Marcott_Next_dT,[colrs(i) '--']);
	%plot(Marcottgl_yr,Marcott_Next_T-Marcott_Next_dT,[colrs(i) '--']);

	%i=3;
	%hpl(i)=plot(HadCRUT_yr,HadCRUT_T,colrs(i)); hold on;

	box off;
	set(hp(2),'yaxisloc','right','color','none');
	set(hp(1),'box','off','xaxisloc','top','xtickl',{''});
	legend(hpl,'GSL','Mann2008','Location','Northwest');
	ylabel('Temp. (C)');
	xlabel('Year CE');
	set(hp,'xlim',[0 2000]);


	pdfwrite(['GSLvstemp' timestr]);
	
	
	% datNAO

	NAOyrs=datNAO.data(:,1);
	NAO=datNAO.data(:,2);

	wl=difftimestep;
	Mavg = abs(bsxfun(@minus,NAOyrs,NAOyrs'))<=wl/2;
	denom=sum(Mavg,2); m=max(denom);
	sub=find(denom==m);
	Mavg=Mavg(sub,:);
	Mavg=bsxfun(@rdivide,Mavg,denom(sub));
	NAOavgyrs = Mavg*NAOyrs;
	NAOavg = Mavg*NAO;

	[intyrs,ia,ib]=intersect(NAOavgyrs,dgradtime);

	k=cov(NAOavg(ia),dgradfs(ib,2));
	scaler = sqrt(k(2,2)/k(1,1))*sign(k(2,2)/k(1,1));
	[c,lags]=xcorr(NAOavg(ia),dgradfs(ib,2));

	lags=-200:200;
	for jj=1:length(lags)
		[intyrs,ia,ib]=intersect(NAOavgyrs,dgradtime-lags(jj));
		clag(jj) = corr(NAOavg(ia),dgradfs(ib,2));
	end
	clf;
%	plot(lags,clag);
%	xlabel('lag');
%	ylabel('correlation');
%	pdfwrite('lagcorr');




	% plot against NAO
	clf;
	subplot(2,1,1);
	hp(1)=subplot(2,1,1)
	colrs='mb';
	clear hpl;
	for j=1:2
		trainsub2 = find((abs(lat-testsites(j,2))<1).*(limiting==0));
		testsub = find(testreg==testsites(j,1));
	%	for i=trainsub2(:)'
	%		plot([time1(i) time2(i)],[1 1]*Y(i)/1000);
	%		hold on;
	%		plot([1 1]*meantime(i),(Y(i)+2*dY0(i)*[-1 1])/1000);
	%	end

		hpl(j)=	plot(testX(testsub,3),1e-3*(f2s(testsub,1)),colrs(j),'linew',2); hold on;
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)+sd2s(testsub,1)),[colrs(j) '--']);
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)-sd2s(testsub,1)),[colrs(j) '--']);
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)+2*sd2s(testsub,1)),[colrs(j) ':']);
		plot(testX(testsub,3),1e-3*(f2s(testsub,1)-2*sd2s(testsub,1)),[colrs(j) ':']);
	end
	xlim([1000 2000]);
	longticks(gca);
	legend(hpl,testnames(1:2),'Location','Southeast');

	ylabel('m');

	clear hpl;
	subplot(2,1,2);
	hpl(1)=plot(dgradtime,dgradfs(:,2),'g','linew',2); hold on;
	plot(dgradtime,dgradfs(:,2)+dgradsds(:,2),'g--');
	plot(dgradtime,dgradfs(:,2)-dgradsds(:,2),'g--');
	plot(dgradtime,dgradfs(:,2)+2*dgradsds(:,2),'g:');
	plot(dgradtime,dgradfs(:,2)-2*dgradsds(:,2),'g:');
	hpl(2)=plot(NAOavgyrs,-scaler*(NAOavg-mean(NAOavg))+mean(dgradfs(:,2)),'r');
	title('NJ - NC (non-linear terms)');
	ylabel(['mm/y (' num2str(difftimestep) '-year avg.)']);
	legend(hpl,'NJ-NC',['-NAO (Trouet) ' num2str(difftimestep) 'y avg'] ,'Location','Southeast');
	xlim([1000 2000]);
	longticks(gca);
	pdfwrite(['NJNCgradvsNAO' timestr]);


end

%

