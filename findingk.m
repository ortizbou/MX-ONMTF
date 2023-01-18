%% Finding the number of common communities, and the number of private and total commmunities per layer.
for r=1:1
    Al=Alr{r};
for l=1:L
    %% Find k_l
    k(l)=EigGap(Al{l});
    %% Find Ul
    [U{l},NMI(l)]=ONMTF(Al{l},GTl{l},k(l));
end

X=[U{:}]';
Z=linkage(X);

kc=0;
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
kp(l)=k(l)-length(find((1+sum(k(1:l-1)))<=Z(1:cut,1:2) & Z(1:cut,1:2)<=sum(k(1:l))));
end

kr{r}=k;
kcr(r)=kc;
kpr{r}=kp;
end