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
## @deftypefn {Function File} {[@var{mesh}]} = MSH3Mstructmesh(@var{x},@var{y},@var{z},@var{region},@var{sides})
##
## Constructs a structured tetrahedral 3D mesh on a parallelepipedal domain,
## and returns a PDEtool-like mesh structure.
##
## Input:
## @itemize @minus
## @item @var{x}: vector representing the 1D meshing of the side
## parallel to x axis.
## @item @var{y}: vector representing the 1D meshing of the side
## parallel to y axis.
## @item @var{z}: vector representing the 1D meshing of the side
## parallel to z axis.
## @item @var{region}: number assigned to the meshed region.
## @item @var{sides}: row vector containing the six numbers assigned to the geometrical surface-edges.
## @end itemize
## 
## Output: mesh basic structure, composed of the following fields
## @itemize @minus
## @item @var{p}: matrix with size 3 times number of mesh point. 
## @itemize @bullet
## @item 1st row: x-coordinates of the points.
## @item 2nd row: y-coordinates of the points.
## @item 3rd row: z-coordinates of the points.
## @end itemize
## @item @var{e}: matrix with size 10 times number of mesh border edges.
## @itemize @bullet
## @item 1st row: p-matrix column number of the first edge-vertex.
## @item 2nd row: p-matrix column number of the second edge-vertex.
## @item 3rd row: p-matrix column number of the third edge-vertex.
## @item 4th row: not initialized, only for compatibility with standard PDE-tool like mesh.
## @item 5th row: not initialized, only for compatibility with standard PDE-tool like mesh.
## @item 6th row: not initialized, only for compatibility with standard PDE-tool like mesh. 
## @item 7th row: not initialized, only for compatibility with standard PDE-tool like mesh. 
## @item 8th row: number of the region to the right of the referred mesh
## edge.
## @item 9th row: number of the region to the left of the referred mesh edge.
## @item 10th row: number of the geometrical border upon which the referred mesh edge is lying on.
## @end itemize
## @item @var{t}:
## @itemize @bullet
## @item 1st row: p-matrix column number of the first tetrahedra vertex.
## @item 2nd row: p-matrix column number of the second tetrahedra vertex.
## @item 3rd row: p-matrix column number of the third tetrahedra vertex.
## @item 4th row: p-matrix column number of the fourth tetrahedra vertex.
## @item 5th row: number of the region upon which the referred trg is lying on.
## @end itemize
## @end itemize 
##
## @seealso{MSH2Mgmsh,MSH2Mjoinstructm,MSH2Msubmesh}
## @end deftypefn

function [mesh] = MSH3Mstructmesh(x,y,z,region,sides)

  ## check for correct input
  scalar = ( isscalar(x) || isscalar(y) || isscalar(z) );
  matrix = ( min(size(x)) + min(size(y)) + min(size(z)) != 3 );
  if scalar
    warning("x, y, z cannot be scalar numbers!")
    print_usage;
  endif
  if matrix
    warning("x, y, z cannot be matrices!")
    print_usage;
  endif

  ## sort point coordinates
  x = sort(x);
  y = sort(y);
  z = sort(z);
  ## compute # of points in each direction
  nx = length(x);
  ny = length(y);
  nz = length(z);

  ## generate verticeces
  [XX,YY,ZZ] = meshgrid(x,y,z);
  p = [XX(:),YY(:),ZZ(:)]';

  iiv (ny,nx,nz)=0;
  iiv(:)=1:nx*ny*nz;
  iiv(end,:,:)=[];
  iiv(:,end,:)=[];
  iiv(:,:,end)=[];
  iiv=iiv(:)';

  ## generate connections:

  ## bottom faces
  n1 = iiv;
  n2 = iiv + 1;
  n3 = iiv + ny;
  n4 = iiv + ny + 1;

  ## top faces
  N1 = iiv + nx * ny;
  N2 = N1  + 1;
  N3 = N1  + ny;
  N4 = N3  + 1;

  t = [...
       [n1; n3; n2; N2],...
       [N1; N2; N3; n3],...
       [N1; N2; n3; n1],...
       [N2; n3; n2; n4],...
       [N3; n3; N2; N4],...
       [N4; n3; N2; n4],...
       ];

  ## generate boundary face list:

  ## left
  T       = t;
  T(:)    = p(1,t)'==x(1);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e1(:,jj)      = t(order(:,ii(jj)),ii(jj));
  endfor
  e1(10,:) = sides(1);

  ## right
  T(:)    = p(1,t)'==x(end);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e2(:,jj)      = t(order(:,ii(jj)),ii(jj));
  end
  e2(10,:) = sides(2);

  ## front
  T(:)    = p(2,t)'==y(1);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e3(:,jj)      = t(order(:,ii(jj)),ii(jj));
  endfor
  e3(10,:) = sides(3);

  ## back
  T(:)    = p(2,t)'==y(end);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e4(:,jj)      = t(order(:,ii(jj)),ii(jj));
  endfor
  e4(10,:) = sides(4);
  
  ## bottom
  T       = t;
  T(:)    = p(3,t)'==z(1);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e5(:,jj)      = t(order(:,ii(jj)),ii(jj));
  endfor
  e5(10,:) = sides(5);
  
  ## top
  T       = t;
  T(:)    = p(3,t)'==z(end);
  [ignore,order] = sort(T,1);
  ii      = (find(sum(T,1)==3));
  order(1,:) = [];
  for jj=1:length(ii)
    e6(:,jj)      = t(order(:,ii(jj)),ii(jj));
  endfor
  e6(10,:) = sides(6);


  mesh.e       = [e1,e2,e3,e4,e5,e6];
  mesh.t       = t;
  mesh.e (9,:) = region;
  mesh.t (5,:) = region;
  mesh.p       = p;

endfunction

%!test
% x = y = z = linspace(0,1,2)
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6)
% assert = (columns(mesh.p),8)
% assert = (columns(mesh.t),6)
% assert = (columns(mesh.e),12)
%!test
% x = y = z = linspace(0,1,3)
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6)
% assert = (columns(mesh.p),27)
% assert = (columns(mesh.t),48)
% assert = (columns(mesh.e),48)
%!test
% x = y = z = linspace(0,1,4)
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6)
% assert = (columns(mesh.p),64)
% assert = (columns(mesh.t),162)
% assert = (columns(mesh.e),108)
%!test
% x = y = z = linspace(0,1,1)
% fail([mesh] = MSH3Mstructmesh(x,y,z,1,1:6),"warning","x, y, z cannot be scalar numbers!")
%!test
% x = y = z = eye(2)
% fail([mesh] = MSH3Mstructmesh(x,y,z,1,1:6),"warning","x, y, z cannot be matrices!")