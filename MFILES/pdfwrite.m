function pdfwrite(figh,fname)

% pdfwrite(figh,fname)
%
% Saves as EPS then converts to PDF using epstopdf.
% (Note that saving as pdf directly puts the figure on a full page.)
%
% INPUTS:
%	figh	figure handle
%	fname	filename (no extension)
%
% Last updated by Bob Kopp rkopp-at-rutgers.edu, 18 May 2014

if nargin==1
        fname=figh;
        figh=gcf;
end

epsfname = [fname '.eps'];
saveas(figh,epsfname,'epsc2');
system(['unset LD_LIBRARY_PATH ; epstopdf ' epsfname]);
delete(epsfname);