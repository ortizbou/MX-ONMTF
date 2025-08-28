%%%%% Multiplex with common communities
%%%%% L layers
function [H_bestNMI,Hl_bestNMI,NMI_realization,averagedNMI,stdNMI, Clusters]=MX_ONMTF(Alr,GTlr,k,kc,kpl)
%  Input: Alr is cell array with r realizations of the networks. Each cell contains another L cells array with the adjacency matrices of each layer
%         k: L vector with the total number of communities per layer
%         kc: number of common communities
%         kpl: L vector with the number of private communities per layer
%         communities membership matrices
%  Output: H_bestNMI and Hl_bestNMI: the embedding matrices of the common and private communities with the best NMI for r realizations
%          NMI_realization, averagedNMI, stdNMI:  NMI values at each realization and the average and standard deviation over these realizations, respectively.
%          Clusters: n*L vector with the labels of the communities
%
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu
%
%
%  Reference: 
%  [1] “Community detection in multiplex networks based on orthogonal nonnegative matrix tri-factorization” Authors: Meiby Ortiz-Bouza and Selin Aviyente

eta=0.5; %% learning rate between (0,1]
realizations=100;

for r=1:realizations
Al=Alr{r};
GTl=GTlr{r};
L=size(Al,2);
n=size(Al{1},2); 
GTSup=vertcat(GTl{:});

%% Running the code multiple times and finding NMI
runs=20;

for j=1:runs
    %% Initializing H1,H2,H,S1,S2,G1,G2
    for l=1:L
    Hl{l}=rand([n,kpl(l)]);
    Sl{l}=diag(rand(kc,1));
    Gl{l}=diag(rand(kpl(l),1));
    end
    H=rand([n,kc]);

    %% Iterative update
    for i=1:1000
        for l=1:L
            Hlnew{l}=Hl{l}.*((Al{l}*Hl{l}*Gl{l} + Hl{l}*Hl{l}'*H*Sl{l}*H'*Hl{l}*Gl{l}')./(H*Sl{l}*H'*Hl{l}*Gl{l} + Hl{l}*Hl{l}'*Al{l}*Hl{l}*Gl{l})).^eta;
            if all(isnan(Hlnew{l}),'all')==1
            Hlnew{l}=rand([n,kpl(l)]);
            end  
        end
        
        numH=0;
        denH=0;
        for l=1:L
        numH=numH + Al{l}*H*Sl{l} + H*H'*Hl{l}*Gl{l}'*Hl{l}'*H*Sl{l};
        denH=denH + Hl{l}*Gl{l}'*Hl{l}'*H*Sl{l} +  H*H'*Al{l}*H*Sl{l};
        end
        Hnew=H.*(numH./denH).^eta;
        if all(isnan(Hnew),'all')==1
        Hnew=rand([n,kc]);
        end

        for l=1:L
            Slnew{l}=Sl{l}.*((H'*Al{l}*H)./(H'*H*Sl{l}*(H'*H) + H'*Hl{l}*Gl{l}*Hl{l}'*H)).^eta;
            if all(isnan(Slnew{l}),'all')==1
            Slnew{l}=diag(rand(kc,1));
            end 

            Glnew{l}=Gl{l}.*((Hl{l}'*Al{l}*Hl{l})./(Hl{l}'*Hl{l}*Gl{l}*(Hl{l}'*Hl{l}) + Hl{l}'*H*Sl{l}*H'*Hl{l})).^eta;
            if all(isnan(Glnew{l}),'all')==1
            Glnew{l}=diag(rand([kpl(l),1]));
            end
        end

        Hres=cellfun(@minus,Hl,Hlnew,'Un',0);
        Sres=cellfun(@minus,Sl,Slnew,'Un',0);
        if (all(norm(Hres{l})<1e-3,'all') && ...
             all(norm(Sres{l})<1e-3,'all') ...
        && all(norm(H-Hnew)<1e-3,'all'))
            break
        end    
    
        Hl=Hlnew;      
        H=Hnew;
        Sl=Slnew;
        Gl=Glnew;
        it=i;
    end

%% Finding Clusters

[Il,ClustersSupra] = assigncomm(Alr,H,Hl,kc,k,kpl,L,'same',1);

%% NMI
for l=1:L
NMIlj{l}=getNMI(GTl{l},Il{l});
end
NMIsupj=getNMI(GTSup,ClustersSupra);
NMI_layerl{j}=NMIlj;
NMI_sup(j)=NMIsupj;
maxNMI=max(NMI_sup);

if maxNMI==NMIsupj
   H_bestNMI{r}=H; 
   Hl_bestNMI{r}=Hl;   
%    Sl_bestNMI{r}=Sl; 
%    Gl_bestNMI{r}=Gl; 
    Clusters=ClustersSupra;
end   

if norm(maxNMI-1)<1e-3
    break
end 


end

[NMI_realization(r),~]= max(NMI_sup);
clear NMI_sup
end

averagedNMI=mean(NMI_realization);

stdNMI=std(NMI_realization);

end

