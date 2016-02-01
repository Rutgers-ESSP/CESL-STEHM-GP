function [wf2,wV2,wsd2,sitespec]=DetrendSLReconstruction(wf,wV,testsites,testreg,testts,firstyears,lastyears,refyear,winwidth)

% [wf2,wV2,wsd2,sitespec]=DetrendSLReconstruction(wf,wV,testsites,testreg,testts,firstyear,lastyear,refyear,winwidth)
%
% Detrend the sea-level reconstruction specified by wf and wV, relative to years
% specified by firstyear and lastyear.
%
% Example:
%
% [GSL,GSLV,GSLsd]=DetrendSLReconstruction(GSL,GSLV,0,testreg,GSL_yrs,0,1800,2000);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 01 14:35:58 EST 2016

%%%

defval('refyear',2000);
defval('firstyears',0);
defval('lastyears',1800);
defval('winwidth',.01);

M=zeros(length(wf),length(wf));
Mavg=eye(length(wf),length(wf));
selfirstyear=ones(size(testsites,1),1)*NaN;
for kk=1:size(testsites,1)
    sub=find((testreg==testsites(kk,1)));
    pp=1;
    sub1=[];
    while (length(sub1)+(pp>length(firstyears)))==0
        if firstyears(pp)~=lastyears(1)
            sub1=find(abs(testts(sub)-firstyears(pp))<=winwidth);
        end
        if length(sub1)==0
            pp=pp+1;
        end 
    end
    sub2=find(abs(testts(sub)-lastyears(1))<=winwidth);
    if (length(sub1)>=1)&&(length(sub2)>=1)
        M(sub,sub(sub1))=-1/length(sub1); M(sub,sub(sub2))=1/length(sub2);
        M(sub,sub)=M(sub,sub)/(mean(testts(sub(sub2)))-mean(testts(sub(sub1))));
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
    Mavg=bsxfun(@eq,testreg,testreg');
    Mavg=bsxfun(@times,Mavg,(testts(:)'>=refyear(1)));    
    Mavg=bsxfun(@times,Mavg,(testts(:)'<=refyear(2)));
    Mavg=bsxfun(@rdivide,Mavg,max(1,sum(Mavg,2)));
    M=(eye(size(M))-Mavg)*M;
end

wf2=M*wf;
wV2=M*wV*M';
wsd2=sqrt(diag(wV2));

subbad=find(ismember(testreg,testsites(find(isnan(selfirstyear)))));
wf2(subbad)=NaN;
wsd2(subbad)=NaN;

sitespec.siteid=testsites;
sitespec.selfirstyear=selfirstyear;
sitespec.trend=trend;