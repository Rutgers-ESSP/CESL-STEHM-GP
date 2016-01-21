function [hl,hk]=PlotWithShadedErrors(dat,rgb,shade,centralline,edgeline,xl,yl)

%  [hl,hk]=PlotWithShadedErrors(dat,rgb,shade,centralline,edgeline,xl,yl)
%
% dat should either be a structure with x, y, dy or a cell array of such structures
%
% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Tue Jul 14 18:51:20 EDT 2015

defval('rgb',[1 0 0 ; 0 1 0 ; 0 0 1 ; 0 0 0]);
defval('shade',.9);
defval('centralline','-');
defval('edgeline','--');
defval('xl',[]); defval('yl',[]);

if isnumeric(shade)
    shades=min(1,rgb+shade);
else
    shades='none';
end

if isstruct(dat)
   dat={dat};
end

for ii=1:length(dat)
    
    
    
    if size(dat{ii}.dy,2)==1
        dy1=dat{ii}.dy;
        dy2=dat{ii}.dy;
    else
        dy1=dat{ii}.dy(:,1);
        dy2=dat{ii}.dy(:,2);
    end
    dat{ii}.dy1=dy1; dat{ii}.dy2=dy2;
    ylow=dat{ii}.y-dy1;
    yhigh=dat{ii}.y+dy2;
    
    xcap=dat{ii}.x;
    
    if length(xl)==2
        xcap=max(xcap,xl(1));
        xcap=min(xcap,xl(2));
    end
    if length(yl)==2
        ylow=max(ylow,yl(1));
        yhigh=min(yhigh,yl(2));
    end

    if isnumeric(shades)
        hk(ii)=patch([xcap ; xcap(end:-1:1)],[yhigh ; ylow(end:-1:1)],shades(ii,:));
        set(hk,'lines','none');
    end
    
    hold on;
end

for ii=1:length(dat)
    if size(dat{ii}.dy,2)==1
        dy1=dat{ii}.dy;
        dy2=dat{ii}.dy;
    else
        dy1=dat{ii}.dy(:,1);
        dy2=dat{ii}.dy(:,2);
    end
    if (length(centralline)>0)&&(~strcmpi(centralline,'none'))
        hl(ii)=plot(dat{ii}.x,dat{ii}.y,centralline,'color',rgb(ii,:),'linew',2); hold on;
    else
        hl=[];
    end
    if (length(edgeline)>0)&&(~strcmpi(edgeline,'none'))
        plot(dat{ii}.x,dat{ii}.y+dat{ii}.dy2,edgeline,'color',rgb(ii,:));
        plot(dat{ii}.x,dat{ii}.y-dat{ii}.dy1,edgeline,'color',rgb(ii,:));
    end
    
    if length(xl)==2
        xlim(xl);
    end
    if length(yl)==2
        ylim(yl);
    end
    
    
end
