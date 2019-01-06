%%  Extract from vectorized FE-Space
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2019-01-06
%
% [varargout] = ffextractfespace(elementNo,VhType,nt,Vh)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  Gets individual FE-Space from a vectorized FE-Space
%

%% Input Parameters
%
%  VhType:    Structure of the FE-Space. For "fespace Vh(Th, [P2,P1])" enter {'P2','P1'}
%  nt:        Number of Meshtriangles
%  vh:        Input FE-Space
%

%% Output Parameters
%
%  varargout: The individual FE-Subspaces
%

%% Licence
%
%  Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
%  This program is free software: you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see
%  <https://www.gnu.org/licenses/>.
%

%% Code

function [varargout] = ffextractfespace(VhType,nt,vh)

    for i=1:numel(VhType)
         %elementType = VhType{elementNo(i)};
         varargout{i} = extract_fespace(i,VhType,nt,vh);
    end

end

%Extract the FE-Space for a particular element from a vectorized FE-Space
function [vhout] = extract_fespace(elementNo,VhType,nt,vhin)

    %DoFs of the implemented vector spaces
    EnDoFs = {{'P1','P2'},{3,6}};
    pattern = [];
    for i=1:numel(VhType)
         idx = strcmpi(VhType{i},EnDoFs{1});
         if any(idx)
             j = find(idx,1);
             %n(i) = EnDoFs{2}{j};
             if i == elementNo
                 pattern = [pattern, true(1,EnDoFs{2}{j})];
             else
                 pattern = [pattern, false(1,EnDoFs{2}{j})];
             end
         else
             fprintf('Element %s is not implemented!\n',VhType{i});
             error('Unknown Element');
         end
    end
    mask = repmat(logical(pattern)',nt,1);
    vhout = vhin(mask);

end

function printhelp()
    fprintf('%s\n\n','Invalid call to ffreadfespace. Correct usage is:');
    fprintf('%s\n',' -- [varargout] = ffreadfespace ()');
end
