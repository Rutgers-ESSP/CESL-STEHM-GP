function varargout=readjxm(fname,degres)
% hmap=READJXM(fname,degres)
% READJXM(fname,degres)
% 
% Reads (or plots, right away) in a spherical harmonics file from Jerry
% Mitrovica and converts it to a full spatial map.
%
% INPUT:
%
% fname        A string with a filename [default: 'seagl_3.92.mn']
% degres       Map resolution [default: Nyquist]
% 
% OUTPUT:
%
% hmap         The map
%
% Last modified by rkopp-at-alumni-dot-caltech-dot-edu, 11/16/2011

% Define defaults
defval('diro','')
defval('fname','seagl_2.1.mn')
defval('degres',[])

% Load the file
fid=fopen(fullfile(diro,fname),'r');
% Read the header
head=fgetl(fid);
LT=sscanf(head,'%f',2);
% Spherical harmonic bandwidth
L=LT(1);
% Time
T=LT(2);
% Read the coefficients
COF=fscanf(fid,'%f');
fclose(fid);
% Check that the length is as I expect it
difer(length(COF)-2*addmup(L)-2*(L/2+1))

% Now rearrange
mods=(L-[0:L]+1);
% The indices of the blank pairs
modu=cumsum(mods+mod(mods,2));
modb=modu(~~mod(mods,2));
blanx=sort([2*modb-1 2*modb]);
% Check these are indeed blanks
difer(COF(blanx))
% Remaining are the non-blanks
COF=skip(COF,blanx);
% Which should be the same number as expected
difer(length(COF)-2*addmup(L))

% Now split into the real and imaginary parts
C=COF(1:2:end);
S=COF(2:2:end);
% The first L+1 should be zero
difer(S(1:L+1))

% Create a blank array with FJS format
[dems,dels,mz,lmcosi]=addmon(L);
modm=[cumsum([0 mods(1:end-1)])+1 addmup(L)+1];
% Populate this with the JXM coefficients
for m=0:L
  % Stick in the "cosine" coefficients
  lmcosi(addmup(m-1:L-1)+m+1,3)=C(modm(m+1):modm(m+2)-1);
  % Stick in the "sine" coefficients
  lmcosi(addmup(m-1:L-1)+m+1,4)=S(modm(m+1):modm(m+2)-1);
end

% Some conventions need to be realigned
CSC=(-1).^dems;
dom=sqrt(2-(dems==0));
lmcosi(:,3)=lmcosi(:,3).*CSC.*dom;
lmcosi(:,4)=lmcosi(:,4).*CSC.*dom;

% The result is not quite right according to our own conventions, but
% we'll take care of that after the expansion rather than adjusting the
% coefficients. The next line could have been done using the
% coefficients, too.
if prod(size(degres))>1
    h=fliplr(plm2xyz(lmcosi,degres(:,1),degres(:,2)));
else
    h=fliplr(plm2xyz(lmcosi,degres));
end
    

if ~nargout
  imagef(h)
  longticks(gca,2)
  plotcont; axis image; kelicol
  caxis([-50 50])
  cb(1)=colorbar('hor');
  longticks(cb,2)
  title(fname)
  fig2print(gcf,'landscape')
  movev(cb,-.1)
  figdisp([],pref(fname))
end

% Provide output
vars={h};
varargout=vars(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=setnan(data)
data(abs(data)<max(abs(data(:)))/100)=NaN;

