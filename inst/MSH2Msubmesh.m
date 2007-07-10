function [omesh,nodelist,elementlist] = MSH2Msubmesh(imesh,intrfc,sdl)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{omesh},@var{nodelist},@var{elementlist}]} = MSH2Msubmesh(@var{imesh},@var{intrfc},@var{sdl})
  ##
  ## Gives as output a specified submesh, and the lists of nodes and element that are mantained.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{imesh}: standard PDEtool-like mesh, with field "p", "e", "t".
  ## @item @var{intrfc}: row vector containing the number of the internal interface sides (numbering referred to mesh.e(5,:)) that should be mantained.
  ##                     Could be empty.
  ## @item @var{sdl} (subdomain list): row vector containing the list of the subdomain that are going to be extracted.
  ## @end itemize
  ##
  ## Output:
  ## @itemize @minus
  ## @item @var{omesh}: standard PDEtool-like mesh, with field "p", "e", "t".
  ## @item @var{nodelist}: list of the node of the original mesh that are present in the restricted one.
  ## @item @var{elementlist}: list of the element of the original mesh that are present in the restricted one.
  ## @end itemize 
  ##
  ## @seealso{MSH2Mstructmesh,MSH2Mjoinstructm,MSH2Mgmsh}
  ## @end deftypefn

  ## This file is part of 
  ##
  ##                   MSH - Meshing Software Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco and Culpo Massimiliano
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
  ##   along with MSH; if not, write to the Free Software
  ##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  ##   USA
  ##
  ##
  ##   MAIN AUTHOR:
  ##   Carlo de Falco
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   AID IN PROGRAMMING AND CLEANING THE CODE: 
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany


  nsd = length(sdl);

  #######
  # set list of output triangles 
  elementlist=[];
  for isd=1:nsd
    elementlist = [elementlist find(imesh.t(4,:) == sdl(isd))];
  end

  omesh.t = imesh.t(:,elementlist);

  #######
  # set list of output nodes
  nodelist          = unique(reshape(imesh.t(1:3,elementlist),1,[]));
  omesh.p         = imesh.p(:,nodelist);


  #######
  # use new node numbering in connectivity matrix
  indx(nodelist) = [1:length(nodelist)];
  iel = [1:length(elementlist)];
  omesh.t(1:3,iel) = indx(omesh.t(1:3,iel));

  #######
  # set list of output edges
  omesh.e =[];
  for isd=1:nsd
    omesh.e = [omesh.e imesh.e(:,imesh.e(7,:)==sdl(isd))];
    omesh.e = [omesh.e imesh.e(:,imesh.e(6,:)==sdl(isd))];
  end
  omesh.e=unique(omesh.e',"rows")';

  #######
  # use new node numbering in boundary segment list
  ied = [1:size(omesh.e,2)];
  omesh.e(1:2,ied) = indx(omesh.e(1:2,ied));

endfunction

%!test
%! [mesh1] = MSH2Mstructmesh(0:.5:1, 0:.5:1, 1, 1:4, 'left');
%! [mesh2] = MSH2Mstructmesh(1:.5:2, 0:.5:1, 1, 1:4, 'left');
%! [mesh] = MSH2Mjoinstructm(mesh1,mesh2,2,4);
%! [omesh,nodelist,elementlist] = MSH2Msubmesh(mesh,[],2);
%! p = [1.00000   1.00000   1.00000   1.50000   1.50000   1.50000   2.00000   2.00000   2.00000
%!      0.00000   0.50000   1.00000   0.00000   0.50000   1.00000   0.00000   0.50000   1.00000];
%! e = [1   1   2   3   4   6   7   8
%!      2   4   3   6   7   9   8   9
%!      0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0
%!      2   5   2   7   5   7   6   6
%!      2   0   2   0   0   0   0   0
%!      1   2   1   2   2   2   2   2];
%! t = [1   2   4   5   2   3   5   6
%!      4   5   7   8   4   5   7   8
%!      2   3   5   6   5   6   8   9
%!      2   2   2   2   2   2   2   2];
%! nl = [7    8    9   10   11   12   13   14   15];
%! el = [9   10   11   12   13   14   15   16];
%! toll = 1e-4;
%! assert(omesh.p,p,toll);
%! assert(omesh.e,e);
%! assert(omesh.t,t);
%! assert(nodelist,nl);
%! assert(elementlist,el);
