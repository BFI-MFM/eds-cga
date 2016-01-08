function formdata(fname,date)
% function formdata(fname,date)
% store endogenous and exogenous variables in a "FRM" TROLL text format file
% INPUT
%   fname: name of the FRM file
%   date:  the date of first observation (i.e. 2007A for an annual dataset)
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none

% Copyright (C) 2007-2010 Dynare Team
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
fid = fopen([fname '_endo.frm'],'w');
n=size(oo_.endo_simul,1);
t=size(oo_.endo_simul,2);
SN=upper(cellstr(M_.endo_names));
for i=1:n
    str=char(SN(i));
    fprintf(fid,'USER: x x DATAFILE: x %s\n',str);
    fprintf(fid,'PER: 1    YEAR: %s   FRAC: 1   NOBS: %d   CLINES: 0   DLINES: ???\n',date,t);
    fprintf(fid,'%10.5f %10.5f %10.5f %10.5f\n',reshape(oo_.endo_simul(i,1:floor(t/4)*4),floor(t/4),4));
    if(floor(t/4)*4<t)
        switch(t-floor(t/4)*4)
          case 1
            fprintf(fid,'%10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
          case 2
            fprintf(fid,'%10.5f %10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
          case 3
            fprintf(fid,'%10.5f %10.5f %10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
        end;
        %else
        %    fprintf(fid,'\n');
    end;
end;
fclose(fid);

fid = fopen([fname '_exo.frm'],'w');
n=size(oo_.exo_simul,2);
t=size(oo_.exo_simul,1);
SN=upper(cellstr(M_.exo_names));
for i=1:n
    str=char(SN(i));
    fprintf(fid,'USER: x x DATAFILE: x %s\n',str);
    fprintf(fid,'PER: 1    YEAR: %s   FRAC: 1   NOBS: %d   CLINES: 0   DLINES: ???\n',date,t);
    fprintf(fid,'%10.5f %10.5f %10.5f %10.5f\n',reshape(oo_.exo_simul(1:floor(t/4)*4,i),floor(t/4),4));
    if(floor(t/4)*4<t)
        switch(t-floor(t/4)*4)
          case 1
            fprintf(fid,'%10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
          case 2
            fprintf(fid,'%10.5f %10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
          case 3
            fprintf(fid,'%10.5f %10.5f %10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
        end;
        %else
        %    fprintf(fid,'\n');
    end;
end;
fclose(fid);
return;