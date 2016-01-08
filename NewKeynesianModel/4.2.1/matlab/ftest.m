function ftest (s1,s2)

% Copyright (C) 2001-2010 Dynare Team
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

global nvx nvy x y lag1

if size(s1,1) ~= 2
    error ('Spécifiez deux fichiers pour la comparaison.') ;
end

for i = 1:2
    if ~ isempty(find(abs(s1(i,:)) == 46))
        error ('Entrez les noms de fichiers sans extensions.') ;
    end
end

s1 = [s1 [' ';' ']] ;
file1 = [s1(1,1:min(find(abs(s1(1,:)) == 32))-1) '.BIN'] ;
file2 = [s1(2,1:min(find(abs(s1(2,:)) == 32))-1) '.BIN'] ;

fid=fopen(file1,'r') ;
n1 = fread(fid,1,'int') ;
n2 = fread(fid,1,'int') ;
n3 = fread(fid,1,'int') ;
lag1 = fread(fid,4,'int') ;
nvx = fread(fid,[n1,n3],'int') ;
x = fread(fid,[n1,n2],'float64') ;
fclose(fid) ;
nvx = char(nvx) ;

fid=fopen(file2,'r') ;
n1 = fread(fid,1,'int') ;
n2 = fread(fid,1,'int') ;
n3 = fread(fid,1,'int') ;
lag2 = fread(fid,4,'int') ;
nvy = fread(fid,[n1,n3],'int') ;
y = fread(fid,[n1,n2],'float64') ;
fclose(fid) ;
nvy = char(nvy) ;

if size(x,1) ~= size(y,1)
    error ('FTEST: The two files don''t have the same number of variables.');
end

for i = 1:size(x,1)
    if ~ strcmp(nvx(i,:),nvy(i,:))
        error ('FTEST: The two files don''t have the same  variables.') ;       
    end
end

if nnz(lag1 - lag2) > 0
    error ('FTEST: Leads and lags aren''t the same in both files.') ;
end

j = zeros(size(s2,1),1);
for i=1:size(s2,1)
    k = strmatch(s2(i,:),nvx,'exact') ;
    if isempty(k)
        t = ['FTEST: Variable ' s2(i) 'doesn''t exist'] ;
        error (t) ;
    else
        j(i) =k;
    end
end

y = y(j,:) ;
x = x(j,:) ;

%06/18/01 MJ replaced beastr by strmatch
