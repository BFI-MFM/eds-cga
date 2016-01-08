function model_info;
%function model_info;

% Copyright (C) 2008-2010 Dynare Team
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

global M_;
fprintf('                                          Informations about %s\n',M_.fname);
fprintf(strcat('                                          ===================',char(ones(1,length(M_.fname))*'='),'\n\n'));
if(isfield(M_,'block_structure'))
    nb_blocks=length(M_.block_structure.block);
    fprintf('The model has %d equations and is decomposed in %d blocks as follow:\n',M_.endo_nbr,nb_blocks);
    fprintf('===============================================================================================================\n');
    fprintf('| %10s | %10s | %30s | %14s | %31s |\n','Block no','Size','Block Type','   Equation','Dependent variable');
    fprintf('|============|============|================================|================|=================================|\n');
    for i=1:nb_blocks
        size_block=length(M_.block_structure.block(i).equation);
        if(i>1)
            fprintf('|------------|------------|--------------------------------|----------------|---------------------------------|\n');
        end;
        for j=1:size_block
            if(j==1)
                fprintf('| %10d | %10d | %30s | %14d | %-6d %24s |\n',i,size_block,Sym_type(M_.block_structure.block(i).Simulation_Type),M_.block_structure.block(i).equation(j),M_.block_structure.block(i).variable(j),M_.endo_names(M_.block_structure.block(i).variable(j),:));
            else
                fprintf('| %10s | %10s | %30s | %14d | %-6d %24s |\n','','','',M_.block_structure.block(i).equation(j),M_.block_structure.block(i).variable(j),M_.endo_names(M_.block_structure.block(i).variable(j),:));
            end;
        end;
    end;
    fprintf('===============================================================================================================\n');
    fprintf('\n');
    for k=1:M_.maximum_endo_lag+M_.maximum_endo_lead+1
        if(k==M_.maximum_endo_lag+1)
            fprintf('%-30s %s','the variable','is used in equations Contemporaneously');
        elseif(k<M_.maximum_endo_lag+1)
            fprintf('%-30s %s %d','the variable','is used in equations with lag ',M_.maximum_endo_lag+1-k);
        else
            fprintf('%-30s %s %d','the variable','is used in equations with lead ',k-(M_.maximum_endo_lag+1));
        end;
        if(size(M_.block_structure.incidence(k).sparse_IM,1)>0)
            IM=sortrows(M_.block_structure.incidence(k).sparse_IM,2);
        else
            IM=[];
        end;
        size_IM=size(IM,1);
        last=99999999;
        for i=1:size_IM
            if(last~=IM(i,2))
                fprintf('\n%-30s',M_.endo_names(IM(i,2),:));
            end;
            fprintf(' %5d',IM(i,1));
            last=IM(i,2);
        end;
        fprintf('\n\n');
    end;
else
    fprintf('There is no block decomposition of the model.\nUse ''block'' model''s option.\n');
end;


function ret=Sym_type(type);
UNKNOWN=0;
EVALUATE_FORWARD=1;
EVALUATE_BACKWARD=2;
SOLVE_FORWARD_SIMPLE=3;
SOLVE_BACKWARD_SIMPLE=4;
SOLVE_TWO_BOUNDARIES_SIMPLE=5;
SOLVE_FORWARD_COMPLETE=6;
SOLVE_BACKWARD_COMPLETE=7;
SOLVE_TWO_BOUNDARIES_COMPLETE=8;
EVALUATE_FORWARD_R=9;
EVALUATE_BACKWARD_R=10;
switch (type)
  case (UNKNOWN),
    ret='UNKNOWN                     ';
  case {EVALUATE_FORWARD,EVALUATE_FORWARD_R},
    ret='EVALUATE FORWARD            ';
  case {EVALUATE_BACKWARD,EVALUATE_BACKWARD_R},
    ret='EVALUATE BACKWARD            ';
  case SOLVE_FORWARD_SIMPLE,
    ret='SOLVE FORWARD SIMPLE        ';
  case SOLVE_BACKWARD_SIMPLE,
    ret='SOLVE BACKWARD SIMPLE        ';
  case SOLVE_TWO_BOUNDARIES_SIMPLE,
    ret='SOLVE TWO BOUNDARIES SIMPLE  ';
  case SOLVE_FORWARD_COMPLETE,
    ret='SOLVE FORWARD COMPLETE      ';
  case SOLVE_BACKWARD_COMPLETE,
    ret='SOLVE BACKWARD COMPLETE      ';
  case SOLVE_TWO_BOUNDARIES_COMPLETE,
    ret='SOLVE TWO BOUNDARIES COMPLETE';
end;




