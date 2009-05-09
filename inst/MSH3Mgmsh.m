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
##  Carlo de Falco <cdf _AT_ users.sourceforge.net>
##
##  Culpo Massimiliano
##  Bergische Universitaet Wuppertal
##  Fachbereich C - Mathematik und Naturwissenschaften
##  Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
##  D-42119 Wuppertal, Germany

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{mesh}]} = MSH3Mgmsh(@var{geometry},@var{option},@var{value},...)
##
## Construct an unstructured 3D mesh making use of the free software gmsh. Give as output the PDE-tool like mesh structure.
##
## Input:
## @itemize @minus
## @item @var{geometry}: name of the ".geo" file describing the 2D geometry. Required by gmsh to start the meshing operation.
## @item @var{option}: option to be used by gmsh
## @item @var{value}: value of the option
## @end itemize
## For more information refer to gmsh manual, or gmsh site:
## 
## http://www.geuz.org/gmsh/
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
## @item 8th row: number of the region to the right of the referred mesh edge.
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

function [mesh] = MSH3Mgmsh(geometry,varargin)

  ## Number of passed options
  noptions  = (nargin - 1) / 2;
  
  optstring = "";
  for ii = 1:noptions
    option    = varargin{2*(ii)-1};
    value     = varargin{2*ii};
    if !ischar(value)
      value = num2str(value);
    endif
    optstring = [optstring," -",option," ",value];
  endfor

  ## Generate mesh using Gmsh
  printf("\n");
  printf("Generating mesh...\n");
  system(["gmsh -format msh -3 -o " geometry ".msh" optstring " " geometry ".geo"]);

  printf("Processing gmsh data...\n");
  ## Points
  com_p   = "awk '/\\$Nodes/,/\\$EndNodes/ {print $2, $3, $4 > ""msh_p.txt""}' ";
  ## Surface edges
  com_e   = "awk '/\\$Elements/,/\\$EndElements/ {if ($2 == ""2"") print $7, $8, $9, $5 > ""msh_e.txt""}' ";
  ## Tetrahedra
  com_t   = "awk '/\\$Elements/,/\\$EndElements/ {if ($2 == ""4"") print $7, $8, $9, $10, $5 > ""msh_t.txt""}' ";
  ## Side edges
  com_s   = "awk '/\\$Elements/,/\\$EndElements/ {if ($2 == ""1"") print $7, $8, $5 > ""msh_s.txt""}' ";

  command = [com_p,geometry,".msh ; "];
  command = [command,com_e,geometry,".msh ; "];
  command = [command,com_t,geometry,".msh ; "];
  command = [command,com_s,geometry,".msh"];
  
  system(command);

  ## Create PDE-tool like structure
  printf("Creating PDE-tool like mesh...\n");
  ## Mesh-points
  p   = load("msh_p.txt")';
  ## Mesh side-edges
  s   = load("msh_s.txt")';
  ## Mesh surface-edges
  tmp = load("msh_e.txt")';
  be  = zeros(10,columns(tmp));
  be([1,2,3,10],:) = tmp;
  ## Mesh tetrahedra
  t   = load("msh_t.txt")';


  ## Remove hanging nodes
  printf("Check for hanging nodes...\n");
  nnodes = columns(p);
  in_msh = intersect( 1:nnodes , t(1:4,:) );
  if length(in_msh) != nnodes
    new_num(in_msh) = [1:length(in_msh)];
    t(1:4,:)        = new_num(t(1:4,:));
    be(1:3,:)       = new_num(be(1:3,:));
    p               = p(:,in_msh);
  endif

  mesh = struct("p",p,"s",s,"e",be,"t",t);

  printf("Deleting temporary files...\n");
  system(["rm -f msh_p.txt msh_e.txt msh_t.txt msh_s.txt *.msh"]);

endfunction