function [I1,ratioI] = matchInt(i1,i0)
S1 = sort(i1(:),'descend');
S0 = sort(i0(:),'descend');
% L0 = round(length(S0)/1e4);
L0 = 50;
ratioI = sum(S0(1:L0))/sum(S1(1:L0));
I1 = i1*ratioI;
end
