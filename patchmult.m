function [pos,q]=patchmult(a,H1,k1)

[~,n]=size(a);
[M1,I1] = max(H1,[],2);

h=zeros(size(H1));

for i=1:k1
h(find(I1==i),i)=1;
end


%%% selecting the column that is not present in layer l
for k=1:k1
[edges,~]=size(find(h(:,k)==1));
di(1)=sum(sum(a.*(h(:,k)*h(:,k)')))/(edges*(edges-1));
[rows,cols]=find((h(:,k)*h(:,k)')==1);
rows=unique(rows);
cols=unique(cols);
m=0;
for i=1:length(rows)
    for j=1:n
        m=m+a(rows(i),j);
    end
end
m=m-sum(sum(a.*(h(:,k)*h(:,k)')));
d0(k)=m/(edges*(n-edges));
end


q=di./d0;


[~,pos]=sort(q);

end
