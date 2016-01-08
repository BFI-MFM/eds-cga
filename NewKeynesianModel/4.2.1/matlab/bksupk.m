function d1 = bksupk(ny,fid,jcf,icc1)

% function d1 = bksupk(ny,fid,jcf,icc1)
% Solves deterministic models recursively by backsubstitution for k leads/lags
%
% INPUTS
%    ny:             number of endogenous variables
%    fid:            saves the elements above the diagonal
%    jcf:            variables index forward
%    icc1:           jacobian column forward
%
% OUTPUTS
%    d1:             vector of backsubstitution results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ options_

icf = [1:jcf-1] ;
ir = [(options_.periods-1)*ny+1:ny*options_.periods] ;
irf = icc1+(options_.periods-1)*ny ;
d1 = zeros(options_.periods*ny,1) ;

ofs = (((options_.periods-1)*ny+1)-1)*jcf*8 ;
junk = fseek(fid,ofs,-1) ;
c = fread(fid,[jcf,ny],'float64')';

d1(ir) = c(:,jcf) ;
ir = ir-ny ;

i = 2 ;

while i <= M_.maximum_lead | i <= options_.periods
    irf1 = selif(irf,irf<=options_.periods*ny) ;

    ofs = (((options_.periods-i)*ny+1)-1)*jcf*8 ;
    junk = fseek(fid,ofs,-1) ;
    c = fread(fid,[jcf,ny],'float64')' ;

    d1(ir) = c(:,jcf) - c(:,1:size(irf1,1))*d1(irf1) ;
    ir = ir - ny ;
    irf = irf - ny ;
    i = i + 1 ;
end

while i <= options_.periods

    ofs = (((options_.periods-i)*ny+1)-1)*jcf*8 ;
    junk = fseek(fid,ofs,-1) ;
    c = fread(fid,[jcf,ny],'float64')' ;

    d1(ir) = c(:,jcf)-c(:,icf)*d1(irf) ;
    ir = ir-ny ;                        
    irf = irf-ny ;
    i = i+1;
end

return ;
