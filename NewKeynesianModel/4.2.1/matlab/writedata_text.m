function writedata_text(fname)
% function writedata(fname)
% store endogenous and exogenous variables in a text file 
% INPUT
%   fname: name of the text file
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none

% Copyright (C) 2007-2009 Dynare Team
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

global M_ oo_
S=[fname '_endo.dat'];
fid = fopen(S,'w');
for i = 1:size(M_.endo_names,1)
    fprintf(fid,'%s ',M_.endo_names(i,:)');
end;
fprintf(fid,'\n');
for i = 1:size(oo_.endo_simul,2)
    fprintf(fid,'%15.7f ',oo_.endo_simul(:,i));
    fprintf(fid,'\n');
end
fclose(fid);

S=[fname '_exo.dat'];
fid = fopen(S,'w');
for i = 1:size(M_.exo_names,1)
    fprintf(fid,'%s ',M_.exo_names(i,:));
end;
fprintf(fid,'\n');
for i = 1:size(oo_.exo_simul,1)
    fprintf(fid,'%15.7f ',oo_.exo_simul(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
return;