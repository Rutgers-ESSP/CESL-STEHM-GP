function A=MergeDataStructures(B,C)

% A=MergeDataStructures(B,C)
%
% Merge data structures B and C.
%
% Last updated by  Bob Kopp, robert-dot-kopp-at-rutgers-dot-edu, Mon Feb 17 23:28:43 EST 2014

fnB=fieldnames(B);
fnC=fieldnames(C);
[fni]=intersect(fnB,fnC);
for ii=1:length(fni)
    if iscell(B.(fni{ii}))
        A.(fni{ii})={B.(fni{ii}){:},C.(fni{ii}){:}};
    elseif size(B.(fni{ii}),1)==size(B.(fni{ii}),2)
        A.(fni{ii})=[B.(fni{ii}) zeros(size(B.(fni{ii}),2),size(C.(fni{ii}),1)) ; zeros(size(C.(fni{ii}),2),size(B.(fni{ii}),1)) C.(fni{ii})];
    else
        A.(fni{ii})=[B.(fni{ii}) ; C.(fni{ii})];
    end
end