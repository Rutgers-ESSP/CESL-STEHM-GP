function [Vnl,Vconst,Vrate,t0,sigmaadj] = PartitionCovarianceLNL(slV,years)

% [Vnl,Vconst,Vrate,t0,sigmaadj] = PartitionCovarianceLNL(slV,years)
% 
% Remove constant rate and constant terms from covariance
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Feb 27 17:06:17 EST 2015
    
    difftimestep=max(diff(years));

    Mdiff = bsxfun(@eq,years,years')-bsxfun(@eq,years,years'+difftimestep);
    sub=find(sum(Mdiff,2)==0);
    Mdiff=Mdiff(sub,:);
    difftimes=bsxfun(@rdivide,abs(Mdiff)*years,sum(abs(Mdiff),2));;
    Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*years);
    
    dV = Mdiff*slV*Mdiff';
    dV=.5*(dV+dV');

    Vrate=max(0,min(dV(:)));
    %    Vnl=slV-Vrate*bsxfun(@times,years-max(years),years'-max(years));
    %Vconst=max(0,min(Vnl(:)));
    %Vnl=Vnl-Vconst;
  
    
    lincv = @(x) x(1)*bsxfun(@times,years-x(2),years'-x(2));
    metric = @(x) sum(sum((slV-lincv([exp(x(1)) x(2)])).^2));

    function [c,ceq] = constr(x)
        delta = (years-x(2)).^2*exp(x(1))-diag(slV);
        c = sum(delta.*(delta>0));
        ceq=0;
    end
 
    t0=2010;
    [x,fval]=fmincon(metric,[log(Vrate) t0],[],[],[],[],[-50 -1000],[8 2010],@constr,optimset('display','off','TolCon',1e-20));
    Vrate=max(1e-10,exp(x(1)));
    Vrate=Vrate.*(Vrate>1e-10);
    
    if constr(x)>0
        [x,fval]=fmincon(metric,[log(Vrate) t0],[],[],[],[],[-50 2010],[8 2010],@constr,optimset('display','off','TolCon',1e-20));
        Vrate=max(1e-10,exp(x(1)));
        Vrate=Vrate.*(Vrate>1e-10);
        
    end        
    
    t0=x(2);
    Vnl = slV-lincv([Vrate t0]);
    sigmaadj=-min([0 ; diag(Vnl)]);
    Vnl=Vnl+eye(size(Vnl))*sigmaadj;
    
    Vconst=max(0,min(Vnl(:)));
    Vnl=Vnl-Vconst;
    
end
