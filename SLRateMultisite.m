function [fslopeavg,Vslopeavg]=SLRateMultisite(wf,wV,testsites,testreg,testts,firstyears,lastyears)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon May 05 15:24:05 EDT 2014

defval('firstyears',[0 1000]);
defval('lastyears',[1800 1800]);
dodiff=((nargout>2)*(length(firstyears)>1));

for pp=1:length(firstyears)
    M=zeros(length(testsites),length(testreg));
    for kk=1:size(testsites,1)
    
        sub=find((testreg==testsites(kk,1)));
        sub1=intersect(sub,find((testts==firstyears(pp))));
        sub2=intersect(sub,find((testts==lastyears(pp))));
        if (length(sub1)==1)&&(length(sub2)==1)
            M(kk,sub1(1))=-1/(lastyears(pp)-firstyears(pp)); M(kk,sub2(1))=1/(lastyears(pp)-firstyears(pp));
            goodsite(kk)=1;
        else
            goodsite(kk)=0;
        end
    end

    fslope=M*wf(:,1);
    Vslope=M*wV(:,:,1)*M';
    sdslope=sqrt(diag(Vslope));

    fslope(find(~goodsite))=NaN;
    sdslope(find(~goodsite))=NaN;
    
    fslopeavg(:,pp) = fslope(:);
    sdslopeavg(:,pp) = sdslope(:);
    Vslopeavg(:,:,pp)=Vslope;


end