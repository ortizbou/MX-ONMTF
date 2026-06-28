function [pos,q]=patchmult(a,H,kc)

%  This function is used in the main algorithm of MX-ONMTF to post-process H and describe how the membership of the common
%  communities is determined across layers, i.e., which columns of H contain information about
%  the common communities present in layer l. This function is explained in Algorithm 3 in [1]
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu
%
%
%  Reference: 
%  [1] “Community detection in multiplex networks based on orthogonal nonnegative matrix tri-factorization” Authors: Meiby Ortiz-Bouza and Selin Aviyente

[~,n]=size(a);
[M1,I1] = max(H,[],2);

h=zeros(size(H));

for i=1:kc
h(find(I1==i),i)=1;
end


%%% selecting the column that is not present in layer l
for k=1:kc
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
