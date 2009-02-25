## Copyright (C) 2007,2008  Carlo de Falco, Massimiliano Culpo
##
## This file is part of 
##
##                   MSH - Meshing Software Package for Octave
## 
##  MSH is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
## 
##  MSH is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with MSH; If not, see <http://www.gnu.org/licenses/>.
##
##
##  AUTHORS:
##  Carlo de Falco
##  Dublin City University
##  School of Mathemetical Sciences
##  Ireland
##
##  Culpo Massimiliano
##  Bergische Universitaet Wuppertal
##  Fachbereich C - Mathematik und Naturwissenschaften
##  Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
##  D-42119 Wuppertal, Germany

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{omesh},@var{nodelist},@var{elementlist}]} = MSH3Msubmesh(@var{imesh},@var{intrfc},@var{sdl})
##
## Gives as output a specified submesh, and the lists of nodes and element that are mantained.
##
## Input:
## @itemize @minus
## @item @var{imesh}: standard PDEtool-like mesh, with field "p", "e", "t".
## @item @var{intrfc}: row vector containing the number of the internal
## interface sides (numbering referred to mesh.e(10,:)) that should be
## mantained. Could be empty.
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
## @seealso{MSH3Mstructmesh,MSH3Mjoinstructm,MSH3Mgmsh}
## @end deftypefn

function [omesh,nodelist,elementlist] = MSH3Msubmesh(imesh,intrfc,sdl)

  ## build element list
  elementlist=[];
  for ir = 1:length(sdl)
    elementlist = [ elementlist find(imesh.t(5,:)==sdl(ir)) ];
  endfor

  ## build nodelist
  nodelist = reshape(imesh.t(1:4,elementlist),1,[]);
  nodelist = unique(nodelist);
  
  ## extract submesh
  omesh.p         = imesh.p  (:,nodelist);
  indx(nodelist)  = 1:length (nodelist);
  omesh.t         = imesh.t  (:,elementlist);
  omesh.t(1:4,:)  = indx(omesh.t(1:4,:));

  omesh.e  = [];
  for ifac = 1:size(imesh.e,2)
    if (length(intersect(imesh.e(1:3,ifac),nodelist) )== 3)
      omesh.e = [omesh.e imesh.e(:,ifac)];
    endif
  endfor

  omesh.e(1:3,:)  = indx(omesh.e(1:3,:));

endfunction

%!shared mesh1,mesh2,jmesh,exmesh,nodelist,elemlist
% x = y = z = linspace(0,1,2);
% x2 = linspace(1,2,2);
% [mesh1] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh2] = MSH3Mstructmesh(x2,y,z,1,1:6);
% [jmesh] = MSH3Mjoinstructm(mesh1,mesh2,2,1);
% [exmesh,nodelist,elemlist] = MSH3Msubmesh(jmesh,2,1);
%!test
% assert(size(exmesh.p),size(mesh1.p))
%!test
% assert(size(exmesh.t),size(mesh1.t))
%!test
% assert(size(exmesh.e),size(mesh1.e))