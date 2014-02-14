% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sat Dec 21 10:03:10 EST 2013



% first optimize only tide gauge data

thetTG2 = thetTGG;
subfixed=subfixedTGG;
ubTG2=ubTGG; lbTG2=lbTGG;

subnotfixed=setdiff(1:length(thetTG2),subfixed);
Mfixed=sparse(length(subnotfixed),length(thetTG2));
for i=1:length(subnotfixed)
	Mfixed(i,subnotfixed(i))=1;
end
fixedvect = 0*thetTG2; fixedvect(subfixed)=thetTG2(subfixed);

trainsub=find(istg==1);
dt1t1=dYears(meantime(trainsub),meantime(trainsub));
dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);
for doglob=[0 1]
	[thetTG2(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x*Mfixed+fixedvect,Ycv0(trainsub,trainsub),dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTG2(subnotfixed),lbTG2(subnotfixed),ubTG2(subnotfixed),doglob);
	disp(sprintf('%0.3f ',thetTG2(subnotfixed)));
end

lbTGG(1:5) = thetTG2(1:5); % set lower bound of amplitudes and temporal scale
thetTGG=thetTG2;

% optimize ignoring geochronological uncertainty

thetTGG2 = [thetTGG .1];
ubTGG2 = [ubTGG 5];
lbTGG2 = [lbTGG 0];

subfixed=union(subfixedTGG,[7 10]); % fix length scales at those determined from the tide gauge data
subnotfixed=setdiff(1:length(thetTGG),subfixed);
Mfixed=sparse(length(subnotfixed),length(thetTGG));
for i=1:length(subnotfixed)
	Mfixed(i,subnotfixed(i))=1;
end
fixedvect = 0*thetTGG; fixedvect(subfixed)=thetTGG2(subfixed);

subnotfixed = [subnotfixed length(thetTGG2)];
trainsub = find((limiting==0));
trainsub=intersect(trainsub,find(meantime>=-1000));
%trainsub=intersect(trainsub,find(datid<1e6));
dt1t1=dYears(meantime(trainsub),meantime(trainsub));
dy1y1 = dDist([lat(trainsub) long(trainsub)],[lat(trainsub) long(trainsub)]);
fp1fp1=bsxfun(@times,obsGISfp(trainsub)-1,obsGISfp(trainsub)'-1);

for doglob=[0]
	[thetTGG2(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv0(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2,dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTGG2(subnotfixed),lbTGG2(subnotfixed),ubTGG2(subnotfixed),doglob);
	disp(sprintf('%0.3f ',thetTGG2(subnotfixed)));
end

% now include geochronological uncertainty, one iteration

wcvfunc = @(x1,x2,thet) cvfuncTGG(x1,x2,dYears(x1,x2),thet,dy1y1,bedmsk(trainsub,trainsub),fp1fp1);

dt = abs(time2-time1)/4;

[dK,df,d2f,yoffset] = GPRdx(meantime(trainsub),Y(trainsub),dt(trainsub),sqrt(dY(trainsub).^2+(thetTGG2(end)*compactcorr(trainsub)).^2),@(x1,x2) wcvfunc(x1,x2,thetTGG),2);

for doglob=[0]
	[thetTGG2(subnotfixed),minima] = SLGPOptimize(Y(trainsub),@(x) traincvTGG(meantime(trainsub),meantime(trainsub),dt1t1,x(1:end-1)*Mfixed+fixedvect,Ycv0(trainsub,trainsub)+diag(x(end)*compactcorr(trainsub)).^2+diag(dK),dy1y1,bedmsk(trainsub,trainsub),fp1fp1),thetTGG2(subnotfixed),lbTGG2(subnotfixed),ubTGG2(subnotfixed),doglob);
	disp(sprintf('%0.3f ',thetTGG2(subnotfixed)));

end

% thetTGG2 = [138.548 37.565 59.492 44.370 360.619 1.195 3.435 0.482 0.178 2.020 0.009 1.000 1.000 0.713 ]
