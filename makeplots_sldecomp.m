% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Dec 21 17:40:33 EST 2013


% SL component plots
% (a) total record
% (b) regional + local record
% (c-d) corresponding average rates of change


for i=1:size(testsites,1)

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

		if (testsites(i,1)==0) && (j==1)
			plot(testX(subA,3),f3s(subA,j)+offsetA,'m');
			plot(testX(subA,3),f3s(subA,j)+sd3s(subA,j)+offsetA,'m--');
			plot(testX(subA,3),f3s(subA,j)-sd3s(subA,j)+offsetA,'m--');
			plot(testX(subA,3),f3s(subA,j)+2*sd3s(subA,j)+offsetA,'m:');
			plot(testX(subA,3),f3s(subA,j)-2*sd3s(subA,j)+offsetA,'m:');
		end
		
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
		xlim([-1000 2010]);
		
%		j=js(1,k);
%		ylim([min(-50,floor(min(slf(subA,j)-2*sqrt(diag(slV(subA,subA,j)))+offsetA)/50)*50) max(50,ceil(max(slf(subA,j)+2*sqrt(diag(slV(subA,subA,j)))+offsetA)/50)*50)]);
		
	end
	
	difftimestep=100;

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
	for n=1:size(noiseMasks,1)
		df2s(:,n)=Mdiff*f2s(:,n);
		dV2s(:,:,n)=Mdiff*V2s(:,:,n)*Mdiff';
		dsd2s(:,n)=sqrt(diag(dV2s(:,:,n)));
	end
	
	clear df3s dV3s dsd3s;
	n=1;
	df3s(:,n)=Mdiff*f3s(:,n);
	dV3s(:,:,n)=Mdiff*V3s(:,:,n)*Mdiff';
	dsd3s(:,n)=sqrt(diag(dV3s(:,:,n)));
	

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

		if (testsites(i,1)==0)&&(j==1)
			plot(difftimes(subA2),df3s(subA2,j),'m');
			plot(difftimes(subA2),df3s(subA2,j)+dsd3s(subA2,j),'m--');
			plot(difftimes(subA2),df3s(subA2,j)-dsd3s(subA2,j),'m--');
			plot(difftimes(subA2),df3s(subA2,j)+2*dsd3s(subA2,j),'m:');
			plot(difftimes(subA2),df3s(subA2,j)-2*dsd3s(subA2,j),'m:');
		end
		
		ylabel(['mm/y (' num2str(difftimestep) '-y)']);
		xlim([-1000 2010]);
		
%		j=js(1,k);
%		ylim([min(-50,floor(min(slf(subA,j)-2*sqrt(diag(slV(subA,subA,j)))+offsetA)/50)*50) max(50,ceil(max(slf(subA,j)+2*sqrt(diag(slV(subA,subA,j)))+offsetA)/50)*50)]);
		
	end
	title(hp(1),testnames2{i})
	longticks(hp);
	[bh,th]=label(hp,'ul',12,[],0,1,1,1.5,1.5); 
	
	pdfwrite(['sldecomp_' testnames{i}])

end

