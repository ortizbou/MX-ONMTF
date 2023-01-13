%%%%% Multiplex with common communities
%%%%% Four layers

% clear all;
% close all;
eta=1; %% learning rate between (0,1]
a=0;
b=0;
c=0;
realizations=1;
[~,n]=size(A{1});
% ModDen1_realization=zeros(1,realizations);
% ModDen2_realization=zeros(1,realizations);
% ModDen3_realization=zeros(1,realizations);
% ModDenSup_realization=zeros(1,realizations);
H_bestModDenSup=cell(6,realizations);
H1_bestModDen=cell(1,realizations);
H1_bestModDenNorm=cell(1,realizations);
H2_bestModDen=cell(1,realizations);
H2_bestModDenNorm=cell(1,realizations);
H3_bestModDen=cell(1,realizations);
H3_bestModDenNorm=cell(1,realizations);
H4_bestModDen=cell(1,realizations);
H4_bestModDenNorm=cell(1,realizations);
H5_bestModDen=cell(1,realizations);
H5_bestModDenNorm=cell(1,realizations);
H6_bestModDen=cell(1,realizations);
H6_bestModDenNorm=cell(1,realizations);

%% Load data

k1=10; %number of communities per layer
k2=10;
k3=10;
k4=10;
k5=10;
K6=10;
kc=10; %number of common communities
kp1=0;
kp2=0;
kp3=0;
kp4=0;
kp5=0;
Kp6=0;
% A1=british;
% A2=airfrance;
%A2=A3;

for r=1:realizations

%% Running the code multiple times and finding NMI
runs=30;

Mod_Den1=-200*ones(1,runs);
Mod_Den2=-200*ones(1,runs);
Mod_Den3=-200*ones(1,runs);
Mod_Den4=-200*ones(1,runs);
Mod_Den5=-200*ones(1,runs);
Mod_Den6=-200*ones(1,runs);
% Mod_Den_Norm1=zeros(1,runs);
% Mod_Den_Norm2=zeros(1,runs);
% Mod_Den_Norm3=zeros(1,runs);
Mod_DenSup=-200*ones(1,runs);

for j=1:runs
    %% Initializing H1,H2,H,S1,S2,G1,G2
    H1=rand([n,kp1]);
    H2=rand([n,kp2]);
    H3=rand([n,kp3]);
    H4=rand([n,kp4]);
    H5=rand([n,kp5]);
    H6=rand([n,kp6]);
    H=rand([n,kc]);
    S1=rand([kc,kc]);   
    S1= tril(S1) + triu(S1', 1);
    S2=rand([kc,kc]);
    S2= tril(S2) + triu(S2', 1);
    S3=rand([kc,kc]);   
    S3= tril(S3) + triu(S3', 1);
    S4=rand([kc,kc]);   
    S4= tril(S4) + triu(S4', 1);
    S5=rand([kc,kc]);   
    S5= tril(S5) + triu(S5', 1);
    S6=rand([kc,kc]);   
    S6= tril(S6) + triu(S6', 1);
    G1=rand([kp1,kp1]);
    G1= tril(G1) + triu(G1', 1);
    G2=rand([kp2,kp2]);
    G2= tril(G2) + triu(G2', 1);
    G3=rand([kp3,kp3]);
    G3= tril(G3) + triu(G3', 1);
    G4=rand([kp4,kp4]);
    G4= tril(G4) + triu(G4', 1);
    G5=rand([kp5,kp5]);
    G5= tril(G5) + triu(G5', 1);
    G6=rand([kp6,kp6]);
    G6= tril(G6) + triu(G6', 1);
    
    
    %% Iterative update
    for i=1:1000
        H1new=H1.*((A1*H1*G1 + H1*H1'*H*S1*H'*H1*G1')./(H*S1*H'*H1*G1 + H1*H1'*A1*H1*G1)).^eta;
%         H1new=H1new/norm(H1new);
%         H1new=[H1new(:,1)/norm(H1new(:,1),1),H1new(:,2)/norm(H1new(:,2),1),H1new(:,3)/norm(H1new(:,3),1)];
        if all(isnan(H1new),'all')==1
        H1new=rand([n,kp1]);
        end  
        H2new=H2.*((A2*H2*G2 + H2*H2'*H*S2*H'*H2*G2')./(H*S2*H'*H2*G2' + H2*H2'*A2*H2*G2)).^eta;  
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
        if all(isnan(H5new),'all')==1
        H5new=rand([n,kp5]);
        end
        H6new=H6.*((A6*H6*G6 + H6*H6'*H*S6*H'*H6*G6')./(H*S6*H'*H6*G6 + H6*H6'*A6*H6*G6)).^eta;
        if all(isnan(H6new),'all')==1
        H6new=rand([n,kp6]);
        end 
        Hnew=H.*((A1*H*S1 + A2*H*S2 + A3*H*S3 +A4*H*S4 + A5*H*S5 + A6*H*S6+ H*H'*H1*G1'*H1'*H*S1 + H*H'*H2*G2'*H2'*H*S2 +H*H'*H3*G3'*H3'*H*S3+ H*H'*H4*G4'*H4'*H*S4+ H*H'*H5*G5'*H5'*H*S5 + H*H'*H6*G6'*H6'*H*S6)./(H1*G1'*H1'*H*S1 ...
            + H2*G2'*H2'*H*S2 +H3*G3'*H3'*H*S3 + H4*G4'*H4'*H*S4 + H5*G5'*H5'*H*S5 + H6*G6'*H6'*H*S6+ H*H'*A1*H*S1 + H*H'*A2*H*S2 + H*H'*A3*H*S3+ H*H'*A4*H*S4+ H*H'*A5*H*S5 + H*H'*A6*H*S6)).^eta;
        
        Hnew=Hnew/norm(Hnew);
        if all(isnan(Hnew),'all')==1
        Hnew=rand([n,kc]);
        end
        S1new=S1.*((H'*A1*H + a*(S2+S3+S4+S5) + b*(1/trace(S1))*eye(kc))./(H'*H*S1*(H'*H) + H'*H1*H1'*H + c*S1)).^eta;
%         S1new=S1new/norm(S1new);
%         S1new=S1new/norm(diag(S1new));
        if all(isnan(S1new),'all')==1
        S1new=diag(rand(kc,1));
        end 
        S2new=S2.*((H'*A2*H + a*(S1+S3+S4+S5+S6) + b*(1/trace(S2))*eye(kc))./(H'*H*S2*(H'*H) + H'*H2*H2'*H + c*S2)).^eta;
%         S2new=S2new/norm(S2new);
%         S2new=S2new/norm(diag(S2new));
        if all(isnan(S2new),'all')==1
        S2new=diag(rand(kc,1)); 
        end 
        S3new=S3.*((H'*A3*H + a*(S2+S1+S4+S5+S6) + b*(1/trace(S3))*eye(kc))./(H'*H*S3*(H'*H) + H'*H3*H3'*H + c*S3)).^eta;
%         S3new=S3new/norm(S3new);
%         S3new=S3new/norm(diag(S3new));
        if all(isnan(S3new),'all')==1
        S3new=diag(rand(kc,1));
        end 
        S4new=S4.*((H'*A4*H + a*(S2+S1+S3+S5+S6) + b*(1/trace(S4))*eye(kc))./(H'*H*S4*(H'*H) + H'*H4*H4'*H + c*S4)).^eta;
        if all(isnan(S4new),'all')==1
        S4new=diag(rand(kc,1));
        end 
        S5new=S5.*((H'*A5*H + a*(S2+S1+S3+S4+S6) + b*(1/trace(S5))*eye(kc))./(H'*H*S5*(H'*H) + H'*H5*H5'*H + c*S5)).^eta;
        if all(isnan(S5new),'all')==1
        S5new=diag(rand(kc,1));
        end 
        S6new=S6.*((H'*A6*H + a*(S2+S1+S3+S4+S5) + b*(1/trace(S6))*eye(kc))./(H'*H*S6*(H'*H) + H'*H6*H6'*H + c*S6)).^eta;
        if all(isnan(S5new),'all')==1
        S6new=diag(rand(kc,1));
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
        if all(isnan(G4new),'all')==1
        G4new=rand([kp4,kp4]);
        end 
        G5new=G5.*((H5'*A5*H5)./(H5'*H5*G5*(H5'*H5) + H5'*H*S5*H'*H5)).^eta;
        if all(isnan(G5new),'all')==1
        G5new=rand([kp5,kp5]);
        end
        G6new=G6.*((H6'*A6*H6)./(H6'*H6*G6*(H6'*H6) + H6'*H*S6*H'*H6)).^eta;
        if all(isnan(G6new),'all')==1
        G6new=rand([kp6,kp6]);
        end 
    
        if (all(norm(H1-H1new)<1e-4,'all') && all(norm(H2-H2new)<1e-4,'all') ...
            && all(norm(H3-H3new)<1e-4,'all') && all(norm(H4-H4new)<1e-4,'all')&&all(norm(H5-H5new)<1e-4,'all') && all(norm(H6-H6new)<1e-4,'all') && all(norm(H-Hnew)<1e-4,'all') && all(norm(S2-S2new)<1e-4,'all') ...
             && all(norm(S3-S3new)<1e-4,'all') && all(norm(S4-S4new)<1e-4,'all') && all(norm(S5-S5new)<1e-4,'all') && all(norm(S6-S6new)<1e-4,'all')  && all(norm(S1-S1new)<1e-4,'all') && all(norm(G1-G1new)<1e-4,'all') ...
            && all(norm(G2-G2new)<1e-4,'all')&& all(norm(G4-G4new)<1e-4,'all')&& all(norm(G5-G5new)<1e-4,'all') && all(norm(G6-G6new)<1e-4,'all'))
            break
        end    
    
        H1=H1new;
        H2=H2new;
        H3=H3new;
        H4=H4new;
        H5=H5new;
        H6=H6new;
        H=Hnew;
        S1=S1new;
        S2=S2new;
        S3=S3new;
        S4=S4new;
        S5=S5new;
        S6=S6new;
        G1=G1new;
        G2=G2new;
        G3=G3new;
        G4=G4new;
        G5=G5new;
        G6=G6new;
        it=i;
    end


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

if H1==inf
    break
end    
if H2==inf
    break
end  
if H3==inf
    break
end  
if H4==inf
    break
end  
if H5==inf
    break
end  
if H6==inf
    break
end 

%% Finding clusters 
% H1=[H,H1];
% H2=[H,H2];
% [M1,I1] = max(H1,[],2);
% [M2,I2] = max(H2,[],2);
% % I2(I2==3)=6;
% % I2(I2==4)=7;
% ClustersSupra=[I1;I2];  

% H1=[H,H1];
% H2=[H,H2];
% H3=[H,H3];
% pos=patchmult(A1,H1,H2,H3,k1);
% H1(:,pos)=0;
% pos=patchmult(A2,H1,H2,H3,k2);
% H2(:,pos)=0;
% pos=patchmult(A3,H1,H2,H3,k3);
% H3(:,pos)=0;

% [M1,I1] = max(H1,[],2);
% [M2,I2] = max(H2,[],2);
% [M2,I3] = max(H3,[],2);

I1=zeros(n,1);
I2=zeros(n,1);
I3=zeros(n,1);
I4=zeros(n,1);
I5=zeros(n,1);
I6=zeros(n,1);
for g=1:kc
    for nodes=1:n
    if all(H(nodes,g)>H1(nodes,:)) && all(H(nodes,g)>H2(nodes,:)) && all(H(nodes,g)>H3(nodes,:)) && all(H(nodes,g)>H4(nodes,:)) && all(H(nodes,g)>H5(nodes,:)) && all(H(nodes,g)>H6(nodes,:))
        I1(nodes)=g;
        I2(nodes)=g;
        I3(nodes)=g;
        I4(nodes)=g;
        I5(nodes)=g;
        I6(nodes)=g;
    end
    end
end    

        nodes=find(I1==0);
        [~,i1(nodes)] = max(H1(nodes,:),[],2);
        I1(nodes)= i1(nodes)+kc;
        nodes=find(I2==0);        
        [~,i2(nodes)] = max(H2(nodes,:),[],2);
        I2(nodes)= i2(nodes)+kc+k1;
        nodes=find(I3==0);        
        [~,i3(nodes)] = max(H3(nodes,:),[],2);
        I3(nodes)= i3(nodes)+kc+k1+k2;
        nodes=find(I4==0);
        [~,i4(nodes)] = max(H4(nodes,:),[],2);
        I4(nodes)= i4(nodes)+kc+k1+k2+k3;
        nodes=find(I5==0);
        [~,i5(nodes)] = max(H5(nodes,:),[],2);
        I5(nodes)= i5(nodes)+kc+k1+k2+k3+k4;
        nodes=find(I6==0);
        [~,i6(nodes)] = max(H6(nodes,:),[],2);
        I6(nodes)= i6(nodes)+kc+k1+k2+k3+k4+k5;


% I1=I1';
% I2=I2';
% I3=I3';

% I2(I2==3)=5;
% I3(I3==3)=6;
% I3(I3==4)=7;
% I3(I3==5)=8;
% I3(I3==6)=9;
% I3(I3==7)=10;
% I4(I4==3)=11;
% I4(I4==4)=12;
% I4(I4==5)=13;
% I5(I5==3)=14;
ClustersSupra=[I1;I2;I3;I4;I5;I6];
H1=[H,H1];
H2=[H,H2];
H3=[H,H3];
H4=[H,H4];
H5=[H,H5];
H6=[H,H6];


%% Modularity density and Modularity Density Normalized


[ModDenj1,ModDenNormj1]=ModDen(kp1+kc,A1,I1);

Mod_Den1(j)=ModDenj1;
maxModDen1=max(Mod_Den1);

if maxModDen1==ModDenj1
   H1_bestModDen{r}=H1; 
end  

% Mod_Den_Norm1(k1,j)=ModDenNormj1;
% maxModDenNorm=max(Mod_Den_Norm1(k1,:));
% 
% if maxModDenNorm==ModDenNormj1
%    H1_bestModDenNorm{r}=H1; 
% end 


[ModDenj2,ModDenNormj2]=ModDen(kp2+kc,A2,I2);

Mod_Den2(j)=ModDenj2;
maxModDen2=max(Mod_Den2);

if maxModDen2==ModDenj2
   H2_bestModDen{r}=H2; 
end  

% Mod_Den_Norm2(k2,j)=ModDenNormj2;
% maxModDenNorm=max(Mod_Den_Norm1(k2,:));
% 
% if maxModDenNorm==ModDenNormj2
%    H2_bestModDenNorm{r}=H2; 
% end 


[ModDenj3,ModDenNormj3]=ModDen(kp3+kc,A3,I3);

Mod_Den3(j)=ModDenj3;
maxModDen3=max(Mod_Den3);

if maxModDen3==ModDenj3
   H3_bestModDen{r}=H3; 
end  

% Mod_Den_Norm3(k3,j)=ModDenNormj3;
% maxModDenNorm=max(Mod_Den_Norm3(k3,:));
% 
% if maxModDenNorm==ModDenNormj3
%    H3_bestModDenNorm{r}=H3; 
% end 

[ModDenj4,ModDenNormj4]=ModDen(kp4+kc,A4,I4);

Mod_Den4(j)=ModDenj4;
maxModDen4=max(Mod_Den4);

if maxModDen4==ModDenj4
   H4_bestModDen{r}=H4; 
end 


[ModDenj5,ModDenNormj5]=ModDen(kp5+kc,A5,I5);

Mod_Den5(j)=ModDenj5;
maxModDen5=max(Mod_Den5);

if maxModDen5==ModDenj5
   H5_bestModDen{r}=H5; 
end  


[ModDenj6,ModDenNormj6]=ModDen(kp6+kc,A6,I6);
Mod_Den6(j)=ModDenj6;
maxModDen6=max(Mod_Den6);

if maxModDen6==ModDenj6
   H6_bestModDen{r}=H6; 
end  
% Mod_DenSup(j)=ModDenj1+ModDenj2+ModDenj3+ModDenj4+ModDenj5;
% maxModDenSup=max(Mod_DenSup);
% 
% if maxModDenSup==ModDenj1+ModDenj2+ModDenj3+ModDenj4+ModDenj5
[ModDenjSup,ModDenNormjSup]=ModDen(kp1+kp2+kp3+kp4+kp5+kp6+kc,[A1 zeros(n,5*n);zeros(n) A2 zeros(n,4*n);zeros(n,2*n) A3 zeros(n,3*n);zeros(n,3*n) A4 zeros(n,2*n); zeros(n,4*n) A5 zeros(n); zeros(n,5*n) A6],ClustersSupra);
Mod_DenSup(j)=ModDenjSup;
maxModDenSup=max(Mod_DenSup);

if maxModDenSup==ModDenjSup
   H_bestModDenSup{1,r}=H1; 
   H_bestModDenSup{2,r}=H2;
   H_bestModDenSup{3,r}=H3;
   H_bestModDenSup{4,r}=H4;
   H_bestModDenSup{5,r}=H5;
   H_bestModDenSup{6,r}=H6;
end  

end

% maxNMI_all(k1)=maxNMI;
% maxModDen1_all(k1)=maxModDen2;
% maxModDen2_all(k1)=maxModDen2;
% maxModDenSup_all(k1)=maxModDenSup;
% maxModDenNorm_all(k1)=maxModDenNorm;

[ModDen1_realization(r),~]=max(Mod_Den1);
[ModDen2_realization(r),~]=max(Mod_Den2);
[ModDen3_realization(r),~]=max(Mod_Den3);
[ModDen4_realization(r),~]=max(Mod_Den4);
[ModDen5_realization(r),~]=max(Mod_Den5);
[ModDen6_realization(r),~]=max(Mod_Den6);

[ModDenSup_realization(r),~]=max(Mod_DenSup);
% [ModDenNorm_realization(r),k1_ModDenNorm(r)]=max(maxModDenNorm_all);

 
end
