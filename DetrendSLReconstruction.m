function [wf2,wV2,wsd2,sitespec]=DetrendSLReconstruction(wf,wV,testsites,testreg,testts,firstyears,lastyears,refyear)

%
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu May 01 23:31:21 EDT 2014

%%%

defval('refyear',2000);
defval('firstyears',0);
defval('lastyears',1800);

M=zeros(length(wf),length(wf));
selfirstyear=ones(length(wf),1)*NaN;
for kk=1:size(testsites,1)
    sub=find((testreg==testsites(kk,1)));
    pp=1;
    sub1=[];
    while (length(sub1)+(pp>length(firstyears)))==0
        sub1=find((testts(sub)==firstyears(pp)));
        if length(sub1)==0
            pp=pp+1;
        end 
    end
    sub2=find((testts(sub)==lastyears(1)));
    if (length(sub1)==1)&&(length(sub2)==1)
        M(sub,sub(sub1(1)))=-1; M(sub,sub(sub2(1)))=1;
        M(sub,sub)=M(sub,sub)/(testts(sub2(1))-testts(sub1(1)));
        selfirstyear(sub) = firstyears(pp);
    end
    trend(kk)=M(sub(1),sub)*wf(sub);
end
M0=M;
M=eye(size(M))-diag(testts-refyear)*M0;

wf2=M*wf;
wV2=M*wV;
wsd2=sqrt(diag(wV2));

subbad=find(isnan(selfirstyear));
wf2(subbad)=NaN;
wsd2(subbad)=NaN;

sitespec.siteid=testsites;
sitespec.selfirstyear=selfirstyear;
sitespec.trend=trend;