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
## @deftypefn {Function File} {[@var{mesh}]} = MSH3Mjoinstructm(@var{mesh1},@var{mesh2},@var{s1},@var{s2})
##
## Join two structured meshes (created by MSH3Mstructmesh) into one
## mesh structure variable.
##
## Input:
## @itemize @minus
## @item @var{mesh1}, @var{mesh2}: standard PDEtool-like mesh, with field "p", "e", "t".
## @item @var{s1}, @var{s2}: number of the corresponding geometrical border edge for respectively mesh1 and mesh2.
## @end itemize
##
## Output:
## @itemize @minus
## @item @var{mesh}: standard PDEtool-like mesh, with field "p", "e", "t".
## @end itemize 
##
## WARNING: on the common edge the two meshes must share the same vertexes.
##
## @seealso{MSH3Mstructmesh,MSH3Mgmsh,MSH3Msubmesh}
## @end deftypefn

function mesh = MSH3Mjoinstructm(mesh1,mesh2,s1,s2)

  ## outside world must always be on the same side of the
  ## boundary of mesh1
  [mesh1.e(8:9,:),I] = sort(mesh1.e(8:9,:));

  ## IF THE REGIONS ARE INVERTED THE VERTEX ORDER SHOULD ALSO BE
  ## INVERTED!!

  ## get interface nodes
  intfcnodes1 = MSH3Mnodesonfaces(mesh1,s1)';
  intfcnodes2 = MSH3Mnodesonfaces(mesh2,s2)';

  ## sort interface nodes by position
  [tmp,I]     = sort(mesh1.p(1,intfcnodes1));
  intfcnodes1 = intfcnodes1(I);
  [tmp,I]     = sort(mesh1.p(2,intfcnodes1));
  intfcnodes1 = intfcnodes1(I);
  [tmp,I]     = sort(mesh1.p(3,intfcnodes1));
  intfcnodes1 = intfcnodes1(I);

  [tmp,I]     = sort(mesh2.p(1,intfcnodes2));
  intfcnodes2 = intfcnodes2(I);
  [tmp,I]     = sort(mesh2.p(2,intfcnodes2));
  intfcnodes2 = intfcnodes2(I);
  [tmp,I]     = sort(mesh2.p(3,intfcnodes2));
  intfcnodes2 = intfcnodes2(I);

  ## delete redundant boundary faces
  ## but first remeber what region
  ## they were connected to
  for is = 1:length(s2)
    ii           = find( mesh2.e(10,:)==s2(is) );
    adreg(is,:)  = unique(mesh2.e(9,ii)); 
  endfor

  for is = 1:length(s2)
    mesh2.e(:,find( mesh2.e(10,:)==s2(is) )) = [];
  endfor

  ## change face numbers
  idx                = [];
  consecutives       = [];
  idx                = unique(mesh2.e(10,:));
  consecutives (idx) = [1:length(idx)] + max(mesh1.e(10,:));
  mesh2.e(10,:)      = consecutives(mesh2.e(10,:));

  ## change node indices in connectivity matrix
  ## and edge list
  idx                   = [];
  consecutives          = [];
  idx                   = 1:size(mesh2.p,2);
  offint                = setdiff(idx,intfcnodes2);
  consecutives (offint) = [1:length(offint)]+size(mesh1.p,2);

  consecutives (intfcnodes2) = intfcnodes1;
  mesh2.e(1:3,:)             = consecutives(mesh2.e(1:3,:));
  mesh2.t(1:4,:)             = consecutives(mesh2.t(1:4,:));

  ## delete redundant points
  mesh2.p(:,intfcnodes2) = [];

  ## set region numbers
  regions             = unique(mesh1.t(5,:));
  newregions(regions) = 1:length(regions);
  mesh1.t(5,:)        = newregions(mesh1.t(5,:));

  ## set region numbers
  regions             = unique(mesh2.t(5,:));
  newregions(regions) = [1:length(regions)]+max(mesh1.t(5,:));
  mesh2.t(5,:)        = newregions(mesh2.t(5,:));

  ## set adjacent region numbers in face structure 2
  [i,j] = find(mesh2.e(8:9,:));
  i    += 7;

  mesh2.e(i,j) = newregions(mesh2.e(i,j));

  ## set adjacent region numbers in edge structure 1
  for is = 1:length(s1)
    ii            = find( mesh1.e(10,:)==s1(is) );
    mesh1.e(8,ii) = newregions(regions(adreg(is,:)));
  endfor

  ## build nbew mesh structure
  mesh.p = [mesh1.p mesh2.p];
  mesh.e = [mesh1.e mesh2.e];
  mesh.t = [mesh1.t mesh2.t];

endfunction

%!shared mesh1,mesh2,jmesh
% x  = y = z = linspace(0,1,2);
% x2 = linspace(1,2,2);
% [mesh1] = MSH3Mstructmesh(x,y,z,1,1:6);
% [mesh2] = MSH3Mstructmesh(x2,y,z,3,1:6);
% [jmesh] = MSH3Mjoinstructm(mesh1,mesh2,2,1);
%!test
% assert(columns(jmesh.p),12)
%!test
% tmp = sort(unique(jmesh.e(10,:)));
% assert(tmp,1:11)
%!test
% assert(columns(jmesh.t),columns(mesh1.t)+columns(mesh2.t))
%!test
% assert(unique(jmesh.e(8:9,:)),0:2)
