function k=EigGap(A)

% realizations=size(Ar,1);
n=size(A,2);
dens=2*sum(A,'all')/(n*(n-1));
rand('seed',100); % reseed so you get a similar picture
G = rand(n,n) < dens;
G = triu(G,1);
G = G + G';
NM=eye(n)-normadj(G);
[~,Dn]=eig(full(NM));
dn=abs(diag(Dn));
dn=sort(dn,'descend');
for i=1:n-1
    dif(i)=dn(i)-dn(i+1);
end    
dif(1)=0;
delta=0.95*max(dif);

% Perform the eigenvalue decomposition
[V,D] = eig(double(A));
gap = zeros(size(A,1),1);
d=sort(abs(diag(D)),'descend');
for i = 1:(size(A,1)-1)
    if d(i) > 0 && d(i+1) > 0
        gap(i) = abs(d(i)) - abs(d(i+1));%/d(i);
    else
        gap(i) = 0;
    end
end
% Choose default number of clusters in the event that maximum of the
% eigengap is non-unique
for j=1:n
if gap(j)>delta
k=j;
end
end

% k=mode(k);
end

