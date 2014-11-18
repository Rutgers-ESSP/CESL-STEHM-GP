function [fslope,sdslope,Vslope,firstyear,lastyear,dosite,fslopediff,sdslopediff,Vslopediff,diffsite,diffplus,diffless]=SLRateCompareCrosssite(wf,wV,testsites,testreg,testX,firstyears,lastyears)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Sun Nov 16 18:07:06 EST 2014

defval('firstyears',[0 1000]);
defval('lastyears',[1800 1800]);
dodiff=((nargout>6)*(length(firstyears)>1));

M=zeros(length(firstyears)*size(testsites,1),length(testreg));
Nrow=0;
for pp=1:length(firstyears)
    for kk=1:size(testsites,1)
        Nrow=Nrow+1;
        sub=find((testreg==testsites(kk,1)));
        sub1=intersect(sub,find(testX(:,3)==firstyears(pp)));
        sub2=intersect(sub,find(testX(:,3)==lastyears(pp)));
        M(Nrow,sub2)=-1/(firstyears(pp)-lastyears(pp));
        M(Nrow,sub1)=1/(firstyears(pp)-lastyears(pp));
        firstyear(Nrow)=firstyears(pp);
        lastyear(Nrow)=lastyears(pp);
        dosite(Nrow)=testsites(kk,1);
    end
end

fslope=M*wf;
Vslope=M*wV*M';
sdslope=sqrt(diag(Vslope));

clear M;
Nrow=0;
if dodiff
    for qq=1:length(firstyears)
        for pp=(qq+1):length(firstyears)
            for kk=1:size(testsites,1)
                Nrow=Nrow+1;
                M(Nrow,:) = zeros(1,length(fslope));
                sub1=find((firstyear==firstyears(qq)).*(lastyear==lastyears(qq)));
                sub2=find((firstyear==firstyears(pp)).*(lastyear==lastyears(pp)));
                M(Nrow,intersect(sub2,find(dosite==testsites(kk,1))))=1;
                M(Nrow,intersect(sub1,find(dosite==testsites(kk,1))))=-1;
                diffsite(Nrow)=testsites(kk,1);
                diffplus(Nrow)=pp;
                diffless(Nrow)=qq;
            end
            
        end
    end
    fslopediff = M*fslope;
    Vslopediff= M*Vslope*M';
    sdslopediff=sqrt(diag(Vslopediff));
    
end
