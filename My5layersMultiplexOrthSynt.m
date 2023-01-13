%%%%% Multiplex with common communities
%%%%% Five layers
function [H_bestNMI,H1_bestNMI,H2_bestNMI,H3_bestNMI,H4_bestNMI,H5_bestNMI,NMI_realization,averagedNMI,stdNMI]=My5layersMultiplexOrthSynt(A1r,GT1r,A2r,GT2r,A3r,GT3r,A4r,GT4r,A5r,GT5r)
tic
% clear all;
% close all;
eta=1; %% learning rate between (0,1]
a=0;
b=0;
c=0;
realizations=2;
% n=256;
%% Load data and GT
% load('Scenario3layers1.mat');      
% load('A1_2.mat');
% A1=A1_block;
% GT1=sort(membership)';
% load('A1_3.mat');
% A2=A1_block;
% GT2=sort(membership)';

% A1=Layer1;
% A2=Layer2;
% H12=X3;

k1=6; %number of communities per layer
k2=6;
k3=5;
k4=4;
k5=4;
kc=3; %number of common communities
kp1=3;
kp2=4;
kp3=3;
kp4=3;
kp5=3;
% GT1=zeros(n,1);
% GT2=zeros(n,1);
% GTSup=zeros(2*n,1);
% GTSup=[GT1;GT2+max(GT1)];
% A1=british;
% A2=airfrance;

for r=1:realizations

A1=full(cell2mat(A1r(r)));
A2=full(cell2mat(A2r(r)));
A3=full(cell2mat(A3r(r)));
A4=full(cell2mat(A4r(r)));
A5=full(cell2mat(A5r(r)));
GT1=full(cell2mat(GT1r(r)));
GT2=full(cell2mat(GT2r(r)));
GT3=full(cell2mat(GT3r(r)));
GT4=full(cell2mat(GT4r(r)));
GT5=full(cell2mat(GT5r(r)));
[~,n]=size(A1);
GTSup=[GT1;GT2;GT3;GT4;GT5];
%% Running the code multiple times and finding NMI
runs=20;
NMI_layer1=zeros(1,runs);
NMI_layer2=zeros(1,runs);
NMI_layer3=zeros(1,runs);
NMI_layer4=zeros(1,runs);
NMI_layer5=zeros(1,runs);
NMI_sup=zeros(1,runs);


for j=1:runs
    %% Initializing H1,H2,H,S1,S2,G1,G2
    H1=rand([n,kp1]);
%     [H1,~] = nndsvd(A1,kp1,0);
%     H1=[H1(:,1)/norm(H1(:,1)),H1(:,2)/norm(H1(:,2)),H1(:,3)/norm(H1(:,3))];
    H2=rand([n,kp2]);
%     [H2,~] = nndsvd(A2,kp2,0);
%      H2=[H2(:,1)/norm(H2(:,1)),H2(:,2)/norm(H2(:,2)),H2(:,3)/norm(H2(:,3)),H2(:,4)/norm(H2(:,4))];
    H3=rand([n,kp3]);
%     [H3,~] = nndsvd(A3,kp3,0);
%      H3=[H3(:,1)/norm(H3(:,1)),H3(:,2)/norm(H3(:,2)),H3(:,3)/norm(H3(:,3))];
    H4=rand([n,kp4]);
%     [H4,~] = nndsvd(A4,kp4,0);
    H5=rand([n,kp5]);
%     [H5,~] = nndsvd(A5,kp5,0);
    H=rand([n,kc]);
%      H=[H(:,1)/norm(H(:,1)),H(:,2)/norm(H(:,2)),H(:,3)/norm(H(:,3))];
    S1=diag(rand(kc,1));
%     S1=eye(kc);
%     S1=rand([kc,kc]);   
%     S1= tril(S1) + triu(S1', 1);
    S2=diag(rand(kc,1));
%     S2=eye(kc);
%     S2=rand([kc,kc]);
%     S2= tril(S2) + triu(S2', 1);
    S3=diag(rand(kc,1));
%     S3=eye(kc);
%     S3=rand([kc,kc]);   
%     S3= tril(S3) + triu(S3', 1);
    S4=diag(rand(kc,1));
%     S4=eye(kc);
    S5=diag(rand(kc,1));
%     S5=eye(kc);
%     G1=eye(kp1);
    G1=diag(rand(kp1,1));
%     G1=rand([kp1,kp1]);
%     G1= tril(G1) + triu(G1', 1);
%     G2=eye(kp2);
    G2=diag(rand(kp2,1));
%     G2=rand([kp2,kp2]);
%     G2= tril(G2) + triu(G2', 1);
%     G3=eye(kp3);
    G3=diag(rand(kp3,1));
%     G3=rand([kp3,kp3]);
%     G3= tril(G3) + triu(G3', 1);
    G4=diag(rand(kp4,1));
%     G4=eye(kp4);
    G5=diag(rand(kp5,1));   
%     G5=eye(kp5);
    
    %% Iterative update
%     tic
    for i=1:1000
        H1new=H1.*((A1*H1*G1 + H1*H1'*H*S1*H'*H1*G1')./(H*S1*H'*H1*G1 + H1*H1'*A1*H1*G1)).^eta;
%         H1new=H1new/norm(H1new);
%         H1new=[H1new(:,1)/norm(H1new(:,1),1),H1new(:,2)/norm(H1new(:,2),1),H1new(:,3)/norm(H1new(:,3),1)];
        if all(isnan(H1new),'all')==1
        H1new=rand([n,kp1]);
        end  
        H2new=H2.*((A2*H2*G2 + H2*H2'*H*S2*H'*H2*G2')./(H*S2*H'*H2*G2 + H2*H2'*A2*H2*G2)).^eta;  
%         H2new=H2new/norm(H2new);
%         H2new=[H2new(:,1)/norm(H2new(:,1),1),H2new(:,2)/norm(H2new(:,2),1),H2new(:,3)/norm(H2new(:,3),1),H2new(:,4)/norm(H2new(:,4),1)];
        if all(isnan(H2new),'all')==1
        H2new=rand([n,kp2]);
        end 
        H3new=H3.*((A3*H3*G3 + H3*H3'*H*S3*H'*H3*G3')./(H*S3*H'*H3*G3 + H3*H3'*A3*H3*G3)).^eta;
%         H3new=H3new/norm(H3new);
%         H3new=[H3new(:,1)/norm(H3new(:,1),1),H3new(:,2)/norm(H3new(:,2),1),H3new(:,3)/norm(H3new(:,3),1)];
        if all(isnan(H3new),'all')==1
        H3new=rand([n,kp3]);
        end  
        H4new=H4.*((A4*H4*G4 + H4*H4'*H*S4*H'*H4*G4')./(H*S4*H'*H4*G4 + H4*H4'*A4*H4*G4)).^eta;
%         H4new=H4new/norm(H4new);
%         H4new=[H4new(:,1)/norm(H4new(:,1),1),H4new(:,2)/norm(H4new(:,2),1),H4new(:,3)/norm(H4new(:,3),1)];
        if all(isnan(H4new),'all')==1
        H4new=rand([n,kp4]);
        end  
        H5new=H5.*((A5*H5*G5 + H5*H5'*H*S5*H'*H5*G5')./(H*S5*H'*H5*G5 + H5*H5'*A5*H5*G5)).^eta;
%         H5new=H5new/norm(H5new);
%         H5new=[H5new(:,1)/norm(H5new(:,1),1),H5new(:,2)/norm(H5new(:,2),1),H5new(:,3)/norm(H5new(:,3),1)];
        if all(isnan(H5new),'all')==1
        H5new=rand([n,kp5]);
        end  
        Hnew=H.*((A1*H*S1 + A2*H*S2 + A3*H*S3 + A4*H*S4 +A5*H*S5 + H*H'*H1*G1'*H1'*H*S1 + H*H'*H2*G2'*H2'*H*S2 +H*H'*H3*G3'*H3'*H*S3 + H*H'*H4*G4'*H4'*H*S4 + H*H'*H5*G5'*H5'*H*S5)./(H1*G1'*H1'*H*S1 ...
            + H2*G2'*H2'*H*S2 +H3*G3'*H3'*H*S3 + H4*G4'*H4'*H*S4 + H5*G5'*H5'*H*S5+ H*H'*A1*H*S1 + H*H'*A2*H*S2 + H*H'*A3*H*S3 + H*H'*A4*H*S4 + H*H'*A5*H*S5)).^eta;
        Hnew=Hnew/norm(Hnew);
        if all(isnan(Hnew),'all')==1
        Hnew=rand([n,kc]);
        end
        S1new=S1.*((H'*A1*H + a*(S2+S3+S4+S5) + b*(1/trace(S1))*eye(kc))./(H'*H*S1*(H'*H) + H'*H1*G1*H1'*H + c*S1)).^eta;
%         S1new=S1new/norm(S1new);
%         S1new=S1new/norm(diag(S1new));
        if all(isnan(S1new),'all')==1
        S1new=diag(rand(kc,1));
        end 
        S2new=S2.*((H'*A2*H + a*(S1+S3+S4+S5) + b*(1/trace(S2))*eye(kc))./(H'*H*S2*(H'*H) + H'*H2*G2*H2'*H + c*S2)).^eta;
%         S2new=S2new/norm(S2new);
%         S2new=S2new/norm(diag(S2new));
        if all(isnan(S2new),'all')==1
        S2new=diag(rand(kc,1)); 
        end 
        S3new=S3.*((H'*A3*H + a*(S2+S1+S4+S5) + b*(1/trace(S3))*eye(kc))./(H'*H*S3*(H'*H) + H'*H3*G3*H3'*H + c*S3)).^eta;
%         S3new=S3new/norm(S3new);
%         S3new=S3new/norm(diag(S3new));
        if all(isnan(S3new),'all')==1
        S3new=diag(rand(kc,1));
        end 
        S4new=S4.*((H'*A4*H + a*(S2+S3+S1+S5)+ b*(1/trace(S4))*eye(kc))./(H'*H*S4*(H'*H) + H'*H4*G4*H4'*H + c*S4)).^eta;
%         S4new=S4new/norm(S4new);
%         S4new=S4new/norm(diag(S4new));
        if all(isnan(S4new),'all')==1
        S4new=diag(rand(kc,1));
        end 
        S5new=S5.*((H'*A5*H + a*(S2+S3+S4+S1)+ b*(1/trace(S5))*eye(kc))./(H'*H*S5*(H'*H) + H'*H5*G5*H5'*H + c*S5)).^eta;
%         S5new=S5new/norm(S5new);
%         S5new=S5new/norm(diag(S5new));
        if all(isnan(S5new),'all')==1
        S5new=diag(rand(kc,1));
        end 
        G1new=G1.*((H1'*A1*H1)./(H1'*H1*G1*(H1'*H1) + H1'*H*S1*H'*H1)).^eta;
%         G1new=G1new/norm(diag(G1new));
        if all(isnan(G1new),'all')==1
        G1new=rand([kp1,kp1]);
        end 
        G2new=G2.*((H2'*A2*H2)./(H2'*H2*G2*(H2'*H2) + H2'*H*S2*H'*H2)).^eta;
%         G2new=G2new/norm(diag(G2new));
        if all(isnan(G2new),'all')==1
        G2new=rand([kp2,kp2]);
        end         
        G3new=G3.*((H3'*A3*H3)./(H3'*H3*G3*(H3'*H3) + H3'*H*S3*H'*H3)).^eta;
%         G3new=G3new/norm(diag(G3new));
        if all(isnan(G3new),'all')==1
        G3new=rand([kp3,kp3]);
        end
        G4new=G4.*((H4'*A4*H4)./(H4'*H4*G4*(H4'*H4) + H4'*H*S4*H'*H4)).^eta;
%         G4new=G4new/norm(diag(G4new));
        if all(isnan(G4new),'all')==1
        G4new=rand([kp4,kp4]);
        end 
        G5new=G5.*((H5'*A5*H5)./(H5'*H5*G5*(H5'*H5) + H5'*H*S5*H'*H5)).^eta;
%         G5new=G5new/norm(diag(G5new));
        if all(isnan(G5new),'all')==1
        G5new=rand([kp5,kp5]);
        end 
    
        if (all(norm(H1-H1new)<1e-4,'all') && all(norm(H2-H2new)<1e-4,'all') ...
            && all(norm(H3-H3new)<1e-4,'all') && all(norm(H4-H4new)<1e-4,'all') && all(norm(H5-H5new)<1e-4,'all') && all(norm(H-Hnew)<1e-4,'all') && all(norm(S2-S2new)<1e-4,'all') ...
             && all(norm(S3-S3new)<1e-4,'all') && all(norm(S4-S4new)<1e-4,'all') && all(norm(S5-S5new)<1e-4,'all') && all(norm(S1-S1new)<1e-4,'all') && all(norm(G1-G1new)<1e-4,'all') ...
           && all(norm(G2-G2new)<1e-4,'all') && all(norm(G4-G4new)<1e-4,'all') && all(norm(G5-G5new)<1e-4,'all'))
            break
        end    
    
        H1=H1new;
        H2=H2new;
        H3=H3new;
        H4=H4new;
        H5=H5new;        
        H=Hnew;
        S1=S1new;
        S2=S2new;
        S3=S3new;
        S4=S4new;
        S5=S5new;
        G1=G1new;
        G2=G2new;
        G3=G3new;
        G4=G4new;
        G5=G5new;
        it=i;
    end

% toc
%% Finding Clusters
% I1=zeros(64,1);
% I2=zeros(64,1);
% for r=1:64
%     if (H1(r,1)<1e-3 && H1(r,2)<1e-3)
%         I1(r)=1;
%     end
% end
% for r=1:64
%     if I1(r)~=1
%         [M1(r),I1(r)] = max(H1(r,:),[],2);
%         I1(r)=I1(r)+1;
%     end    
% end    
% 
% for r=1:64
%     if (H2(r,1)<1e-3)
%         I2(r)=1;
%     end
% end
% for r=1:64
%     if I2(r)~=1
%         [M2(r),I2(r)] = max(H2(r,:),[],2);
%         I2(r)=I2(r)+3;
%     end    
% end   
% imagesc(H)
% if H1==inf
%     break
% end    
% if H2==inf
%     break
% end    
% if H3==inf
%     break
% end 
% if H4==inf
%     break
% end 
% if H5==inf
%     break
% end 


H1=[H,H1];
H2=[H,H2];
H3=[H,H3];
H4=[H,H4];
H5=[H,H5];

pos=patchmult(A2,H2,k2);
H2(:,pos(1))=0;

pos=patchmult(A3,H3,k3);
H3(:,pos(1))=0;

pos=patchmult(A4,H4,k4);
H4(:,pos(1))=0;
H4(:,pos(2))=0;

pos=patchmult(A5,H5,k5);
H5(:,pos(1))=0;
H5(:,pos(2))=0;

[~,I1] = max(H1,[],2);
[~,I2] = max(H2,[],2);
[~,I3] = max(H3,[],2);
[~,I4] = max(H4,[],2);
[~,I5] = max(H5,[],2);
% I2(I2==3)=8;
I2(I2==4)=10;
I2(I2==5)=11;
I2(I2==6)=12;
I2(I2==7)=13;
% I3(I3==3)=12;
I3(I3==4)=16;
I3(I3==5)=17;
I3(I3==6)=18;
% I4(I4==3)=16;
I4(I4==4)=21;
I4(I4==5)=22;
I4(I4==6)=23;
% I5(I5==3)=19;
I5(I5==4)=25;
I5(I5==5)=26;
I5(I5==6)=27;
ClustersSupra=[I1;I2;I3;I4;I5]; 

% for g=1:kc
%     for nodes=1:n
%     if all(H(nodes,g)>H1(nodes,:)) && all(H(nodes,g)>H2(nodes,:))
%         I1(nodes)=g;
%         I2(nodes)=g+k1;
%     else 
%         [~,i1(nodes)] = max(H1(nodes,:),[],2);
%         I1(nodes)= i1(nodes)+kc;
%         [~,i2(nodes)] = max(H2(nodes,:),[],2);
%         I2(nodes)= i2(nodes)+kc+k1;
%     end 
%     end
% end
% 
% ClustersSupra=[I1';I2'];
% H1=[H,H1];
% H2=[H,H2];

%% NMI

NMI1j=getNMI(GT1,I1);
NMI2j=getNMI(GT2,I2);
NMI3j=getNMI(GT3,I3);
NMI4j=getNMI(GT4,I4);
NMI5j=getNMI(GT5,I5);
NMIsupj=getNMI(GTSup,ClustersSupra);
NMI_layer1(j)=NMI1j;
NMI_layer2(j)=NMI2j;
NMI_layer3(j)=NMI3j;
NMI_layer4(j)=NMI4j;
NMI_layer5(j)=NMI5j;
maxNMI1=max(NMI_layer1);
maxNMI2=max(NMI_layer2);
maxNMI3=max(NMI_layer3);
maxNMI4=max(NMI_layer4);
maxNMI5=max(NMI_layer5);
NMI_sup(j)=NMIsupj;
maxNMI=max(NMI_sup);


if maxNMI==NMIsupj
   H_bestNMI{r}=H;  
   H1_bestNMI{r}=H1; 
   H2_bestNMI{r}=H2;
   H3_bestNMI{r}=H3;
   H4_bestNMI{r}=H4;
   H5_bestNMI{r}=H5;   
%    S1_bestNMI{r}=S1; 
%    S2_bestNMI{r}=S2;
%    S3_bestNMI{r}=S3;
%    S4_bestNMI{r}=S4;
%    S5_bestNMI{r}=S5;
%    G1_bestNMI{r}=G1; 
%    G2_bestNMI{r}=G2;
%    G3_bestNMI{r}=G3;
%    G4_bestNMI{r}=G4;  
%    G5_bestNMI{r}=G5;
end   

if maxNMI==1
    break
end 

if maxNMI==1.000
    break
end 


end



[NMI_realization(r),k1_NMI(r)]= max(NMI_sup);

end

averagedNMI=mean(NMI_realization);


stdNMI=std(NMI_realization);

toc

end


