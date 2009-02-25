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
## @deftypefn {Function File} {[@var{nodelist}]} =@
## MSH3Mnodesonfaces(@var{mesh}, @var{facelist})
##
## Returns a list of the nodes lying on the faces @var{facelist} of
## the mesh @var{mesh}.
##
## Input:
## @itemize @minus
## @item @var{mesh}: standard PDEtool-like mesh, with field "p", "e", "t".
## @item @var{facelist}: row vector containing the number of the faces (numbering referred to mesh.e(10,:)).
## @end itemize
##
## Output:
## @itemize @minus
## @item @var{nodelist}: list of the nodes that lies on the specified faces.
## @end itemize 
##
## @seealso{MSH3Mgeomprop,MSH3Mtopprop}
## @end deftypefn

function [nodelist] = MSH3Mnodesonfaces(mesh,facelist);

  facefaces = [];
  
  for ii=1:length(facelist)
    facefaces = [facefaces,find(mesh.e(10,:)==facelist(ii))];
  endfor

  facenodes = mesh.e(1:3,facefaces);
  nodelist  = unique(facenodes(:));
  
endfunction

%!shared x,y,z,mesh
% x = y = z = linspace(0,1,2);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6);
%!test
% nodelist = MSH3Mnodesonfaces(mesh,1);
% assert(nodelist,[1 2 5 6]')
%!test
% nodelist = MSH3Mnodesonfaces(mesh,2);
% assert(nodelist,[3 4 7 8]')
%!test
% nodelist = MSH3Mnodesonfaces(mesh,3);
% assert(nodelist,[1 3 5 7]')
%!test
% nodelist = MSH3Mnodesonfaces(mesh,[1 2 3]);
% assert(nodelist,[1:8]')