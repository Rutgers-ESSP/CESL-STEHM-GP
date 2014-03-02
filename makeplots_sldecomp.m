function makeplots_sldecomp(dataset,f2s,sd2s,V2s,testloc,labl)

% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Feb 18 23:13:46 EST 2014



% SL component plots
% (a) total record
% (b) regional + local record
% (c-d) corresponding average rates of change

defval('labl','');
xlim0=[-1000 2010];

angd= @(Lat0,Long0,lat,long) (180/pi)*(atan2(sqrt((cosd(lat).*sind(long-Long0)).^2+(cosd(Lat0).*sind(lat)-sind(Lat0).*cosd(lat).*cosd(long-Long0)).^2),(sind(Lat0).*sind(lat)+cosd(Lat0).*cosd(lat).*cosd(long-Long0))));

dYears=@(years1,years2) abs(bsxfun(@minus,years1',years2));
dDist=@(x1,x2)angd(repmat(x1(:,1),1,size(x2,1)),repmat(x1(:,2),1,size(x2,1)),repmat(x2(:,1)',size(x1,1),1),repmat(x2(:,2)',size(x1,1),1))'+1e6*(bsxfun(@plus,x1(:,1)',x2(:,1))>1000);

testreg=testloc.reg;
testsites=testloc.sites;
testnames=testloc.names;
testnames2=testloc.names2;
testX=testloc.X;

istg = dataset.istg;
meantime=dataset.meantime;
lat=dataset.lat;
long=dataset.long;
Y=dataset.Y0;
Ycv0=dataset.Ycv0;
limiting=dataset.limiting;
obsGISfp=dataset.obsGISfp;
compactcorr=dataset.compactcorr;
time1=dataset.time1;
time2=dataset.time2;
dY=dataset.dY;

for i=1:size(testsites,1)
    wxlim=[];

	figure;
	subA = find(testreg == testsites(i,1));
	distfrom=dDist(testsites(i,2:3),[lat long]);
	subB=find(distfrom<.1);
	
	clf; clear hp;
	
	if testsites(i,1)>0
		js = [1 3 4];
	else
		js=[1 6];
	end
	
	for k=1:length(js)
		offsetA=0;
		j=js(1,k);
		hp(k)=subplot(2,length(js),k);
		box on;
		
		hold on
		plot(testX(subA,3),f2s(subA,j)+offsetA,'r');
		plot(testX(subA,3),f2s(subA,j)+sd2s(subA,j)+offsetA,'r--');
		plot(testX(subA,3),f2s(subA,j)-sd2s(subA,j)+offsetA,'r--');
		plot(testX(subA,3),f2s(subA,j)+2*sd2s(subA,j)+offsetA,'r:');
		plot(testX(subA,3),f2s(subA,j)-2*sd2s(subA,j)+offsetA,'r:');

		if k==1
			for nn=1:length(subB)
				plot([1 1]*meantime(subB(nn)),Y(subB(nn))+[-1 1]*dY(subB(nn)));
				plot([time1(subB(nn)) time2(subB(nn))],Y(subB(nn))*[1 1]);
			end
		end
		
		if j==3
			title('Regional + Local');
		elseif j==4
			title('Regional + Local non-linear');
		elseif j==6
			title('Greenland');
		end


		ylabel('mm');
		if length(wxlim)=0
    		wxlim = xlim0;
    		wxlim(1) = max(floor(min(testX(subA,3))/100)*100,xlim0(1));
        end
		xlim(wxlim);
				
	end
	
	difftimestep=100;

	Mdiff = bsxfun(@eq,testX(:,3),testX(:,3)')-bsxfun(@eq,testX(:,3),testX(:,3)'+difftimestep);
	Mdiff = Mdiff .* bsxfun(@eq,testreg,testreg');
	sub=find(sum(Mdiff,2)==0);
	Mdiff=Mdiff(sub,:);
	difftimes=bsxfun(@rdivide,abs(Mdiff)*testX(:,3),sum(abs(Mdiff),2));;
	diffreg=bsxfun(@rdivide,abs(Mdiff)*testreg,sum(abs(Mdiff),2));;
	Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testX(:,3));

	clear df2s dV2s dsd2s;
	for n=1:size(f2s,2)
		df2s(:,n)=Mdiff*f2s(:,n);
		dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
		dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
	end
	
	subA2 = find(diffreg == testsites(i,1));

	for k=1:length(js)
		offsetA=0;
		hp(k+length(js))=subplot(2,length(js),k+length(js));
		box on;
		
		hold on
		j=js(1,k);
		plot(difftimes(subA2),df2s(subA2,j),'r');
		plot(difftimes(subA2),df2s(subA2,j)+dsd2s(subA2,j),'r--');
		plot(difftimes(subA2),df2s(subA2,j)-dsd2s(subA2,j),'r--');
		plot(difftimes(subA2),df2s(subA2,j)+2*dsd2s(subA2,j),'r:');
		plot(difftimes(subA2),df2s(subA2,j)-2*dsd2s(subA2,j),'r:');

		ylabel(['mm/y (' num2str(difftimestep) '-y)']);
		xlim(wxlim);
		
	end
	title(hp(1),testnames2{i})
	longticks(hp);
	[bh,th]=label(hp,'ul',12,[],0,1,1,1.5,1.5); 
	
	pdfwrite(['sldecomp_' testnames{i} labl])

end

