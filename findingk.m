%% Finding the number of common communities, and the number of private and total commmunities per layer.
L=size(Al,2);

for l=1:L
    %% Find k_l
    k(l)=EigGap(Al{l});
    %% Find Ul
%     [U{l},NMI(l)]=ONMTF(Al{l},GTl{l},k(l)); %%% for synthetic networks
    U{l}=ONMTF_real(Al{l},k(l)); %%% for real networks
end

X=[U{:}]';
Z=linkage(X);

kc=0;
cut=size(Z,1);
for i=2:size(X,1)-1 
    if max(Z(i-1,1:2))<=size(X,1)
        kc=kc+1;
    end
    d(i)=(Z(i,3)-Z(i-1,3))/Z(i-1,3);
    if d(i)>=0.5
    cut=i;
    break
    end
end

for l=1:L
kpl(l)=k(l)-length(find((1+sum(k(1:l-1)))<=Z(1:cut,1:2) & Z(1:cut,1:2)<=sum(k(1:l))));
end
