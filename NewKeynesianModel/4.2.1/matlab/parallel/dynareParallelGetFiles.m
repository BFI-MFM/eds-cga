function dynareParallelGetFiles(NamFileInput,PRCDir,Parallel)
% PARALLEL CONTEXT
% In a parallel context, this is a specialized mono-directional (Remote to Local) version of copy()
% function.
%
%
% INPUTS
%  o NamFileInput   []   ...
%  o PRCDir         []   ...
%  o Parallel       []   ...
%
%  OUTPUTS
%  None
%
%
%
% Copyright (C) 2009-2011 Dynare Team
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


if ischar(NamFileInput),
    for j=1:size(NamFileInput,1),
        NamFile(j,:)={['.',filesep],deblank(NamFileInput(j,:))};
    end
    NamFileInput = NamFile;
end

for indPC=1:length(Parallel),
    if Parallel(indPC).Local==0,
        if ~ispc, %isunix || (~matlab_ver_less_than('7.4') && ismac),
            for jfil=1:size(NamFileInput,1),
                % if ~isempty(dynareParallelDir(NamFileInput{jfil,2},[PRCDir,filesep,NamFileInput{jfil,1}],Parallel(indPC))),

                if exist('OCTAVE_VERSION') % Patch for peculiar behaviour of ls under Linux.
                    % It is necessary to manage the jolly char '*'!

                    FindAst=strfind(NamFileInput{jfil,2},'comp_status_random_walk_metropolis_hastings_core*');
                    
                    if isempty (FindAst)

                        [NonServeL NonServeR]= system(['scp ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,':',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',NamFileInput{jfil,1}]);

                    else

                        filenameTemp=NamFileInput{jfil,2};
                        
                        [NotUsed FlI]=system(['ssh ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,' ls ',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',filenameTemp, ' 2> OctaveStandardOutputMessage.txt']);
                        
                        if isempty (FlI)
                            return
                        end

                        AstPos=strfind(filenameTemp,'.mat')-1;
                        FiMat=findstr(FlI, '.mat');
                        NumFileToCopy=length(FiMat);


                        for i=1: NumFileToCopy
                            Ni=num2str(i);
                            filenameTemp(1,AstPos)=Ni;
                            [NonServeL NonServeR]= system(['scp ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,':',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1},filenameTemp,' ',NamFileInput{jfil,1}]);
                        end
                    end

                else

                    [NonServeL NonServeR]= system(['scp ',Parallel(indPC).UserName,'@',Parallel(indPC).ComputerName,':',Parallel(indPC).RemoteDirectory,'/',PRCDir,'/',NamFileInput{jfil,1},NamFileInput{jfil,2},' ',NamFileInput{jfil,1}]);
                end


                % end
            end
        else
            for jfil=1:size(NamFileInput,1),
                if ~isempty(dynareParallelDir(NamFileInput{jfil,2},[PRCDir,filesep,NamFileInput{jfil,1}],Parallel(indPC))),
                    copyfile(['\\',Parallel(indPC).ComputerName,'\',Parallel(indPC).RemoteDrive,'$\',Parallel(indPC).RemoteDirectory,'\',PRCDir,'\',NamFileInput{jfil,1},NamFileInput{jfil,2}],NamFileInput{jfil,1});
                end
            end
        end
    end
end
