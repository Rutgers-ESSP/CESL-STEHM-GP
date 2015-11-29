function [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(wf,wV,testsites,testreg,testts,firstyears,lastyears)

% [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]=SLRateCompare(wf,wV,testsites,testreg,testts,firstyears,lastyears)
%
% Calculate sea-level rates and differences in rates within sites, over the time intervals
% defined by firstyears and lastyears.
%
% INPUTS:
%
%     wf - column vector of sea-level heights
%     wV - a corresponding covariance matrix
%     testsites - a matrix of site information, with the first column being site id
%     testreg - column vector listing the siteid for each element of wf
%     testts - column vector of ages for each element of wf
%     firstyears - row vector of interval starts
%     lastyears - row vector of interval ends
%
% OUTPUTS:
%
%     fslopeavg - matrix of rates, with rows corresponding to sites
%                 and columns to intervals
%     sdslopeavg - corresponding matrix of standard deviations
%     fslopeavgdiff - matrix of interperiod rate differences, with
%                     rows corresponding to sites
%     sdslopeavgdiff - corresponding standard deviations
%     diffplus - the interval being added to calculate differences
%     diffless - the intervla being substracted to calculate differences
%
% EXAMPLE:
%
% firstyears=[0 700 1400];
% lastyears= [700 1400 1800];
%
% [fslopeavg,sdslopeavg,fslopeavgdiff,sdslopeavgdiff,diffplus,diffless]= ...
%      SLRateCompare(f2s{1}(:,1),V2s{1}(:,:,1),testsites,testreg, ...
%      testX(:,3),firstyears,lastyears);
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Wed Nov 05 22:41:45 EST 2014

defval('firstyears',[0 1000]);
defval('lastyears',[1800 1800]);
dodiff=((nargout>2)*(length(firstyears)>1));

for kk=1:size(testsites,1)
    
    sub=find((testreg==testsites(kk,1)));
    M=zeros(length(firstyears),length(sub));
    for pp=1:length(firstyears)

        sub1=find((testts(sub)==firstyears(pp)));
        sub2=find((testts(sub)==lastyears(pp)));
        if (length(sub1)==1)&&(length(sub2)==1)
            M(pp,sub1(1))=-1; M(pp,sub2(1))=1;
            M(pp,:)=M(pp,:)/(lastyears(pp)-firstyears(pp));
        end
    end

    fslope=M*wf(sub,1);
    Vslope=M*wV(sub,sub,1)*M';
    sdslope=sqrt(diag(Vslope));

    notgood=find(sum(abs(M),2)==0);
    modifier=0*sdslope; modifier(notgood)=1e6;
    sdslope=sdslope+modifier; Vslope=Vslope+diag(modifier.^2);
    
    if dodiff
        counter=1;
        for pp=1:length(firstyears)
            for qq=(pp+1):length(firstyears)
                Mdiff(counter,pp)=-1;
                Mdiff(counter,qq)=1;
                diffplus(counter)=qq; diffless(counter)=pp;
                counter=counter+1;
            end
        end
        fslopediff = Mdiff*fslope;
        Vslopediff = Mdiff*Vslope*Mdiff';
        sdslopediff = sqrt(diag(Vslopediff));

        fslopediff(find(sdslopediff>1e4))=NaN;
        sdslopediff(find(sdslopediff>1e4))=NaN;
        
        fslopeavgdiff(kk,:) = fslopediff(:)';
        sdslopeavgdiff(kk,:)= sdslopediff(:)';
    end
    
    fslope(find(sdslope>1e4))=NaN;
    sdslope(find(sdslope>1e4))=NaN;
    
    fslopeavg(kk,:) = fslope(:)';
    sdslopeavg(kk,:) = sdslope(:)';


end

fslope(find(sdslope>1e4))=NaN;
sdslope(find(sdslope>1e4))=NaN;

if dodiff
    fslopediff(find(sdslopediff>1e4))=NaN;
    sdslopediff(find(sdslopediff>1e4))=NaN;
else
    fslopeavgdiff=[];
    sdslopeavgdiff=[];
    diffplus=[];
    diffless=[];
end