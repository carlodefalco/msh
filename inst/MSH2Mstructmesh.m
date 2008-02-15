function [mesh] = MSH2Mstructmesh(x,y,region,sides,varargin)
  
  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{mesh}]} = MSH2Mstructmesh(@var{x},@var{y},@var{region},@var{sides},@var{string})
  ##
  ## Constructs a structured 2D mesh on a rectangular domain, and gives as output the PDE-tool mesh structure.
  ##
  ## Input:
  ## @itemize @minus
  ## @item @var{x}: vector representing the 1D meshing of the side parallel to x axis.
  ## @item @var{y}: vector representing the 1D meshing of the side parallel to y axis.
  ## @item @var{region}: number assigned to the meshed region.
  ## @item @var{sides}: row vector containing the four numbers assigned to the geometrical edges.
  ## @item @var{string}: (optional) orientation of the diagonal edge of the structured mesh.
  ##                     Values: "right", "left", "random". Default is "right".
  ## @end itemize
  ##
  ## Output: mesh basic structure, composed of the following fields
  ## @itemize @minus
  ## @item @var{p}: matrix with size 2 times number of mesh point. 
  ## @itemize @bullet
  ## @item 1st row: x-coordinates of the points.
  ## @item 2nd row: y-coordinates of the points.
  ## @end itemize
  ## @item @var{e}: matrix with size 7 times number of mesh border edges.
  ## @itemize @bullet
  ## @item 1st row: p-matrix column number of the first edge-vertex.
  ## @item 2nd row: p-matrix column number of the second edge-vertex.
  ## @item 3rd row: not initialized, only for compatibility with standard PDE-tool like mesh.
  ## @item 4th row: not initialized, only for compatibility with standard PDE-tool like mesh.
  ## @item 5th row: number of the geometrical border upon which the referred mesh edge is lying on.
  ## @item 6th row: number of the region to the right of the referred mesh edge.
  ## @item 7th row: number of the region to the left of the referred mesh edge.
  ## @end itemize
  ## @item @var{t}:
  ## @itemize @bullet
  ## @item 1st row: p-matrix column number of the first trg-vertex.
  ## @item 2nd row: p-matrix column number of the second trg-vertex.
  ## @item 3rd row: p-matrix column number of the third trg-vertex.
  ## @item 4th row: number of the region upon which the referred trg is lying on.
  ## @end itemize
  ## @end itemize 
  ##
  ## @seealso{MSH2Mgmsh,MSH2Mjoinstructm,MSH2Msubmesh}
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
  ##   along with MSH; If not, see <http://www.gnu.org/licenses/>.
  ##
  ##
  ##   MAIN AUTHORS:
  ##   Carlo de Falco
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany

  default = 'right';

  if length(varargin)==0
    string = default;
  else
    string = varargin{1};
  end

  switch string
  case 'right'
    [mesh] = Ustructmesh_right(x, y, region, sides);
  case 'left'
    [mesh] = Ustructmesh_left(x, y, region, sides);
  case 'random'
    [mesh] = Ustructmesh_random(x, y, region, sides);
  otherwise
    error(['Passed string has not a valid value: ' string]);
  endswitch

endfunction

##RIGHT DIAGONAL STRUCTURED MESH
function [mesh]=Ustructmesh_right(x,y,region,sides)
  
  x = sort(x);
  y = sort(y);

  nx = length(x);
  ny = length(y);
  [XX,YY] = meshgrid(x,y);
  p = [XX(:),YY(:)]';
  iiv (ny,nx)=0;
  iiv(:)=1:nx*ny;
  iiv(end,:)=[];
  iiv(:,end)=[];
  iiv=iiv(:)';
  t = [[iiv;iiv+ny;iiv+ny+1],[iiv;iiv+ny+1;iiv+1] ];
  t (4,:)=region;

  l1 = 1+ny*([1:nx]-1);
  l4 = 1:ny;
  l2 = ny*(nx-1)+1:nx*ny;
  l3 = ny + l1 -1;

  e = [ l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])
       l1([2:end]) l2([2:end]) l3([2:end]) l4([2:end])
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       l1([1:end-1])*0+sides(1) l2([1:end-1])*0+sides(2) l3([1:end-1])*0+sides(3) l4([1:end-1])*0+sides(4)
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0+region
	   ];
  mesh.p = p; mesh.e = e; mesh.t = t;

endfunction

##LEFT DIAGONAL STRUCTURED MESH
function [mesh]=Ustructmesh_left(x,y,region,sides)

  x = sort(x);
  y = sort(y);

  nx = length(x);
  ny = length(y);
  [XX,YY] = meshgrid(x,y);
  p = [XX(:),YY(:)]';
  iiv (ny,nx)=0;
  iiv(:)=1:nx*ny;
  iiv(end,:)=[];
  iiv(:,end)=[];
  iiv=iiv(:)';
  t = [[iiv;iiv+ny;iiv+1],[iiv+1;iiv+ny;iiv+ny+1] ];
  t (4,:)=region;

  l1 = 1+ny*([1:nx]-1);
  l4 = 1:ny;
  l2 = ny*(nx-1)+1:nx*ny;
  l3 = ny + l1 -1;

  e = [ l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])
       l1([2:end]) l2([2:end]) l3([2:end]) l4([2:end])
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       l1([1:end-1])*0+sides(1) l2([1:end-1])*0+sides(2) l3([1:end-1])*0+sides(3) l4([1:end-1])*0+sides(4)
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0+region
	   ];
  mesh.p = p; mesh.e = e; mesh.t = t;

endfunction

##RANDOM DIAGONAL STRUCTURED MESH
function [mesh]=Ustructmesh_random(x,y,region,sides)

  x = sort(x);
  y = sort(y);

  nx = length(x);
  ny = length(y);
  [XX,YY] = meshgrid(x,y);
  p = [XX(:),YY(:)]';
  iiv (ny,nx)=0;
  iiv(:)=1:nx*ny;
  iiv(end,:)=[];
  iiv(:,end)=[];
  iiv=iiv(:)';

  niiv = length(iiv);
  theperm = iiv(randperm(niiv));
  first = theperm(1:floor(niiv/2));
  second = theperm(floor(niiv/2)+1:end);

  t = [[first;first+ny;first+ny+1],[first;first+ny+1;first+1] ];
  t = [t,[second;second+ny;second+1],[second+ny;second+ny+1;second+1] ];

  t (4,:)=region;

  l1 = 1+ny*([1:nx]-1);
  l4 = 1:ny;
  l2 = ny*(nx-1)+1:nx*ny;
  l3 = ny + l1 -1;

  e = [ l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])
       l1([2:end]) l2([2:end]) l3([2:end]) l4([2:end])
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
       l1([1:end-1])*0+sides(1) l2([1:end-1])*0+sides(2) l3([1:end-1])*0+sides(3) l4([1:end-1])*0+sides(4)
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0
	   [l1([1:end-1]) l2([1:end-1]) l3([1:end-1]) l4([1:end-1])]*0+region
	   ];
  mesh.p = p; mesh.e = e; mesh.t = t;

endfunction
