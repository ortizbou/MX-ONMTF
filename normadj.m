function A = normadj(A)
%% compute the normalized adjacency matrix

% This code is extracted from the file sgwt_laplacian.m, which is part of
% the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

N=size(A,1);
degrees=vec(full(sum(A)));
% to deal with loops, must extract diagonal part of A
diagw=diag(A);

% w will consist of non-diagonal entries only
[ni2,nj2,w2]=find(A);
ndind=find(ni2~=nj2); % as assured here
ni=ni2(ndind);
nj=nj2(ndind);
w=w2(ndind);

di=vec(1:N); % diagonal indices  

% normalized laplacian D^(-1/2)*(D-A)*D^(-1/2)
% diagonal entries
dL=(diagw./degrees); % will produce NaN for degrees==0 locations
dL(degrees==0)=0;% which will be fixed here
% nondiagonal entries
ndL=w./vec( sqrt(degrees(ni).*degrees(nj)) );
L=sparse([ni;di],[nj;di],[ndL;dL],N,N);
A=full(L);
