function [pos,q]=patchmult(a,H1,k1)

[~,n]=size(a);
[M1,I1] = max(H1,[],2);

h=zeros(size(H1));

for i=1:k1
h(find(I1==i),i)=1;
end


%%% selecting the column that is not present in layer 2
[edges,~]=size(find(h(:,1)==1));
di(1)=sum(sum(a.*(h(:,1)*h(:,1)')))/(edges*(edges-1));
[rows,cols]=find((h(:,1)*h(:,1)')==1);
rows=unique(rows);
cols=unique(cols);
% patch0=ones(size(h1(:,1)*h1(:,1)'))-h1(:,1)*h1(:,1)';
m=0;
for i=1:length(rows)
    for j=1:n
        m=m+a(rows(i),j);
    end
end
m=m-sum(sum(a.*(h(:,1)*h(:,1)')));
d0(1)=m/(edges*(n-edges));

[edges,~]=size(find(h(:,2)==1));
di(2)=sum(sum(a.*(h(:,2)*h(:,2)')))/(edges*(edges-1));
[rows,cols]=find((h(:,2)*h(:,2)')==1);
rows=unique(rows);
cols=unique(cols);
m=0;
for i=1:length(rows)
    for j=1:n
        m=m+a(rows(i),j);
    end
end 
m=m-sum(sum(a.*(h(:,2)*h(:,2)')));
d0(2)=m/(edges*(n-edges));

[edges,~]=size(find(h(:,3)==1));
di(3)=sum(sum(a.*(h(:,3)*h(:,3)')))/(edges*(edges-1));
[rows,cols]=find((h(:,3)*h(:,3)')==1);
rows=unique(rows);
cols=unique(cols);
m=0;
for i=1:length(rows)
    for j=1:n
        m=m+a(rows(i),j);
    end
end    
m=m-sum(sum(a.*(h(:,3)*h(:,3)')));
d0(3)=m/(edges*(n-edges));

q=di./d0;


[~,pos]=sort(q);

% [edges,~]=size(find(h1(:,2)==1));
% s1(1)=sum(sum(A1.*(h1(:,2)*h1(:,2)')))/(edges*(edges-1));
% [edges,~]=size(find(h1(:,2)==1));
% s1(2)=sum(sum(A2.*(h1(:,2)*h1(:,2)')))/(edges*(edges-1));
% [edges,~]=size(find(h1(:,2)==1));
% s1(3)=sum(sum(A3.*(h1(:,2)*h1(:,2)')))/(edges*(edges-1));


 end
