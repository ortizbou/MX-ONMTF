%% Applying Symmetric Nonnegative Matrix 3-Factorization with orthogonal restriction to each layer separately
%%%% 50 runs are performed and the U1 with the best NMI is chosen. This is
%%%% repeated over 100 realizations of the same network and then averaged
%%%% Modularity is reported. Try different values of k

function [U1_bestMod]=ONMTF_real(A,k)


%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu


runs=50;
eta=0.5; 

% A=Alr{1};
[~,n]=size(A);

%% Running the code multiple times and finding NMI

for j=1:runs
    %% Initializing U1,U2
U1=rand([n,k]);
C1=diag(rand([k,1]));

%% Iterative update
for i=1:1000
    U1new=U1.*((A*U1*C1')./(U1*U1'*A*U1*C1')).^eta; 
    C1new=C1.*((U1'*A*U1)./(U1'*U1*C1*(U1'*U1))).^eta;
    if all(isnan(U1new),'all')==1
        U1new=rand([n,k1]);
    end  
    if all(isnan(C1new),'all')==1
        C1new=rand([k1,k1]);
    end  
    if (all(norm(U1-U1new)<1e-6,'all') && all(norm(C1-C1new)<1e-6,'all'))
        break
    end  
    U1=U1new;
    C1=C1new;    
end
    
%% Finding Clusters and Computing Modularity Density
[~,I1] = max(U1,[],2);
ModDenj=ModDen(k,A,I1);
ModD(j)=ModDenj;
maxModD=max(ModD);


if maxModD==ModDenj 
   U1_bestMod=U1; 
end   

 
end 


end
