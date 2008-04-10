function [nodelist] = MSH2Mnodesonsides(mesh,sidelist)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{nodelist}]} =@
  ## MSH2Mnodesonsides(@var{mesh}, @var{sidelist})
  ##
  ## Returns a list of the nodes lying on the sides @var{sidelist} of
  ## the mesh @var{mesh}.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{mesh}: standard PDEtool-like mesh, with field "p", "e", "t".
  ## @item @var{sidelist}: row vector containing the number of the sides (numbering referred to mesh.e(5,:)).
  ## @end itemize
  ##
  ## Output:
  ## @itemize @minus
  ## @item @var{nodelist}: list of the nodes that lies on the specified sides.
  ## @end itemize 
  ##
  ## @seealso{MSH2Mgeomprop,MSH2Mtopprop}
  ## @end deftypefn

  ## This file is part of 
  ##
  ##                   MSH - Meshing Software Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco
  ##              Copyright (C) 2007  Culpo Massimiliano
  ## 
  ##   MSH is free software; you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation; either version 2 of the License, or
  ##   (at your option) any later version.
  ## 
  ##   MSH is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ## 
  ##   You should have received a copy of the GNU General Public License
  ##   along with MSH; If not, see <http://www.gnu.org/licenses/>.
  ##
  ##
  ##   AUTHORS:
  ##   Carlo de Falco
  ##   Dublin City University
  ##   School of Mathemetical Sciences
  ##   Ireland
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universitaet Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
  ##   D-42119 Wuppertal, Germany



  edgelist    =[];

  for ii = 1:length(sidelist)
    edgelist=[edgelist,find(mesh.e(5,:)==sidelist(ii))];
  end

  ##Set list of nodes with Dirichelet BCs
  nodelist = mesh.e(1:2,edgelist);
  nodelist = [nodelist(1,:) nodelist(2,:)];
  nodelist = unique(nodelist);

endfunction

%!test
%! [mesh1] = MSH2Mstructmesh(0:.5:1, 0:.5:1, 1, 1:4, 'left');
%! [mesh2] = MSH2Mstructmesh(1:.5:2, 0:.5:1, 1, 1:4, 'left');
%! [mesh] = MSH2Mjoinstructm(mesh1,mesh2,2,4);
%! [nodelist] = MSH2Mnodesonsides(mesh,[1 2]);
%! reallist = [1   4   7   8   9];
%! assert(nodelist,reallist);
