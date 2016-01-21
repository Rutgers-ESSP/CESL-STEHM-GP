function [fdiff,sddiff,Vdiff,diffpairs,diffplus,diffless,difft]=SLSiteCompare(wf,wV,testsites,testreg,testX)

% [fdiff,sddiff,Vdiff,diffpairs,diffplus,diffless,difft]=SLSiteCompare(wf,wV,testsites,testreg,testX)
%
% Compare differences in sea levels between sites.
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Thu Jul 16 10:45:11 EDT 2015

counter=0;
diffplus=[];
diffless=[];
difft=[];
diffpairs=[];
testts=testX(:,3);
for kk=1:(size(testsites,1)-1)
    for jj=(kk+1):size(testsites,1)
        counter=counter+1;
        sub1=find(testreg==testsites(kk));
        sub2=find(testreg==testsites(jj));
        [u,ui,uj]=intersect(testts(sub1),testts(sub2));
        M=sparse(length(u),length(wf));
        M(:,sub1(ui))=-1;
        M(:,sub2(uj))=1;
        M=M.*bsxfun(@eq,u,testts');
        if (counter>1)
            Mfull=[Mfull ; M];
        else
            Mfull=M;
        end
        diffplus=[diffplus ; repmat(testsites(jj),size(M,1),1)];
        diffless=[diffless ; repmat(testsites(kk),size(M,1),1)];
        difft=[difft ; u(:)];
        diffpairs=[diffpairs ; testsites(jj) testsites(kk)];
    end
end

fdiff=Mfull*wf;
Vdiff=Mfull*wV*Mfull';
sddiff=sqrt(diag(Vdiff));

