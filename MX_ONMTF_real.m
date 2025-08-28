%%%%% Multiplex with common communities
%%%%% L layers
function [H_best, Hl_best, Clusters]=MX_ONMTF_real(Alr,k,kc,kpl)
%  Input: Alr is cell array with r realizations of the networks. Each cell contains another L cells array with the adjacency matrices of each layer
%         k: L vector with the total number of communities per layer
%         kc: number of common communities
%         kpl: L vector with the number of private communities per layer
%         communities membership matrices
%  Output: H_best is the R realization H with the best metric at each realization
%          Clusters: n*L vector with the labels of the communities
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu

eta=0.5; %% learning rate between (0,1]
realizations=1;

for r=1:realizations
Al=Alr;
L=size(Al,2);
n=size(Al{1},2); 

%% Running the code multiple times and finding NMI
runs=20;

for j=1:runs
    %% Initializing Hl,H,Sl,Gl
    for l=1:L
    Hl{l}=rand([n,kpl(l)]);
    Sl{l}=diag(rand(kc,1));
    Gl{l}=diag(rand(kpl(l),1));
    end
    H=rand([n,kc]);

    %% Iterative update
    for i=1:2000
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

[Il,ClustersSupra] = assigncomm(Alr,H,Hl,kc,k,kpl,L,'same',0.8);  %% choose the threshold for the assignment of teh common communities


% %% Modularity
% for l=1:L
% [ModDenj{l},~]=ModDen(k(l),Al{l},Il{l});
% end
% Asup=zeros(L*n);
% for l=1:L
%     Asupra(1+n*(l-1):n*l,1+n*(l-1):n*l)=Al{l};
% end
% ModDenSupj=ModDen(sum(kpl)+kc,Asupra,ClustersSupra);
% Mod_DenSup(j)=ModDenSupj;
% maxModDenSup=max(Mod_DenSup);
% 
% if maxModDenSup==ModDenSupj
%    H_best{r}=H; 
%    Hl_best{r}=Hl;   
%    Sl_best{r}=Sl; 
% %    Gl_best{r}=Gl; 
%    Clusters=ClustersSupra;
% end
%% Trace minimization
for l=1:L
    tr(l)=norm(Al{l}-(H*Sl{l}*H'+Hl{l}*Gl{l}*Hl{l}'),'fro')+trace(Hl{l}'*Hl{l}-eye(kpl(l)));
end
trall(j)=sum(tr)+trace(H'*H-eye(kc));
mintrace=min(trall);

if mintrace==trall(j)  %%%% change j to r and remove %
   H_best{r}=H; 
   Hl_best{r}=Hl;   
   Sl_best{r}=Sl; 
%    Gl_best{r}=Gl; 
   Clusters{r}=ClustersSupra;
end

end

% [ModDen_realization(r),~]= max(Mod_DenSup);
% clear Mod_DenSup

% clear trall tr %%% change if using more realizations
end


end

