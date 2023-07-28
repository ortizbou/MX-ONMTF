function [ModularityDensity,ModularityDensityNorm]=ModDen(k,A,I)
 
%%% Input k: number of communities, I: cluster labels, A: Adjacency matrix
% This function computes the multilayer Modularity Density and Normalized Modularity Density of a community partition of the network using the labels and the adjacency matrix
% Modularity Density is used to quatify a community structure when the ground truth (true labels) are unknown.
% 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------

 
%  LV=zeros(1,k);
%  LVn=zeros(1,k);
ModularityDensity=0;
ModularityDensityNorm=0;
ModularityDensity_new=0;
ModularityDensityNorm_new=0;
absV=0;
absVn=0;
% row=0;
% rown=0;
LV=0;
LVn=0;
LV_new=0;
LVn_new=0;

for c=1:k
    LV=0;
    LVn=0;
    LV_new=0;
    LVn_new=0;
row = find(I==c);
if isempty(row)
   ModularityDensity_new=ModularityDensity_new;
   ModularityDensityNorm_new=ModularityDensityNorm_new;
else
 absV=length(row);
 for l= 1:length(row)
     for h= 1:length(row)
         LV_new= LV+ A(row(l),row(h));
         LV=LV_new;
     end
 end
 
  rown = find(I~=c);
  if isempty(rown)
      LVn=0;
      absVn=1;
  else    
  absVn=length(rown);
    for l= 1:length(row)
        for h=1:length(rown)
         LVn_new= LVn+ A(row(l),rown(h));
         LVn=LVn_new;
        end
    end
  end

 ModularityDensity_new=ModularityDensity + (LV-LVn)/absV;
 ModularityDensityNorm_new=ModularityDensityNorm + (LV/(absV*(absV-1)) - LVn/(absV*absVn));
 ModularityDensity=ModularityDensity_new;
 ModularityDensityNorm=ModularityDensityNorm_new;
end
end
 
% ModularityDensitySum=sum(ModularityDensity);
% ModularityDensityNormSum=sum(ModularityDensityNorm);

end
