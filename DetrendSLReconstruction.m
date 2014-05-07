function [wf2,wV2,wsd2,sitespec]=DetrendSLReconstruction(wf,wV,testsites,testreg,testts,firstyears,lastyears,refyear)

%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 05 08:32:28 EDT 2014

%%%

defval('refyear',2000);
defval('firstyears',0);
defval('lastyears',1800);

M=zeros(length(wf),length(wf));
selfirstyear=ones(size(testsites,1),1)*NaN;
for kk=1:size(testsites,1)
    sub=find((testreg==testsites(kk,1)));
    pp=1;
    sub1=[];
    while (length(sub1)+(pp>length(firstyears)))==0
        if firstyears(pp)~=lastyears(1)
            sub1=find((testts(sub)==firstyears(pp)));
        end
        if length(sub1)==0
            pp=pp+1;
        end 
    end
    sub2=find((testts(sub)==lastyears(1)));
    if (length(sub1)==1)&&(length(sub2)==1)
        M(sub,sub(sub1(1)))=-1; M(sub,sub(sub2(1)))=1;
        M(sub,sub)=M(sub,sub)/(testts(sub2(1))-testts(sub1(1)));
        selfirstyear(kk) = firstyears(pp);
    end
    trend(kk)=M(sub(1),sub)*wf(sub);
    if isnan(selfirstyear(kk))
        trend(kk)=NaN;
    end
end

M0=M;
M=eye(size(M))-diag(testts-refyear(end))*M0;
if length(refyear)==2
    Mavg = (testts(:)'>=firstyears).*(testts(:)'<=firstyears);
    Mavg=Mavg/sum(Mavg);
    Mavg=repmat(Mavg,length(testts),1);
    M=(eye(size(M))-Mavg)*M;
end

wf2=M*wf;
wV2=M*wV;
wsd2=sqrt(diag(wV2));

subbad=find(ismember(testreg,testsites(find(isnan(selfirstyear)))));
wf2(subbad)=NaN;
wsd2(subbad)=NaN;

sitespec.siteid=testsites;
sitespec.selfirstyear=selfirstyear;
sitespec.trend=trend;