## Copyright (C) 2007,2009  Carlo de Falco, Massimiliano Culpo
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
## @deftypefn {Function File} {[@var{mesh}]} = MSH2Mgmsh(@var{geometry},@var{option},@var{value},...)
##
## Construct an unstructured triangular 2D mesh making use of the free
## software gmsh. Return a PDE-tool like mesh structure.
##
## Input:
## @itemize @minus
## @item @var{geometry}: basename of the ".geo" file to be meshed.
## @item @var{option}: string of the option to pass gmsh.
## @item @var{value}: value of the option to pass gmsh.
## @end itemize
##
## For more information regarding the possible option to pass, refer to gmsh manual or gmsh site:
## http://www.geuz.org/gmsh/
##
## @var{mesh} structure is composed of the following fields:
##
## @itemize @minus
## @item @var{p}: 2 X (# nodes) matrix. Contain mesh points coordinates. 
## @item @var{e}: 7 X (# side edges) matrix. Contain mesh side
## edges information.
## @item @var{t}: 4 X (# triangles) matrix. Contain pointer to @var{p}
## field, as well as region number.
## @end itemize 
##
## @seealso{MSH2Mstructmesh, MSH2Mjoinstructm, MSH2Msubmesh}
## @end deftypefn

function [mesh] = MSH2Mgmsh(geometry,varargin)

  ## 1. Check input
  ## 1a. Number of input
  if !mod(nargin,2)
    warning("WRONG NUMBER OF INPUT.");
    print_usage;
  endif

  ## 2. Build mesh
  noptions  = (nargin - 1) / 2; # Number of passed options
  
  ## 2a. Construct system command string
  optstring = "";
  for ii = 1:noptions
    option    = varargin{2*(ii)-1};
    value     = varargin{2*ii};
    if !ischar(value)
      value = num2str(value);
    endif
    optstring = [optstring," -",option," ",value];
  endfor

  ## 2b. Invoke gmsh
  printf("\n");
  printf("Generating mesh...\n");
  system(["gmsh -format msh -2 -o " geometry ".msh" optstring " " geometry ".geo"]);

  ## 2c. Build structure fields
  printf("Processing gmsh data...\n");
  ## Points
  com_p   = "awk '/\\$Nodes/,/\\$EndNodes/ {print $2, $3 > ""msh_p.txt""}' ";
  ## Side edges
  com_e   = "awk '/\\$Elements/,/\\$EndElements/ {if ($2 == ""1"") print $7, $8,$5 > ""msh_e.txt""}' ";
  ## Triangles
  com_t   = "awk '/\\$Elements/,/\\$EndElements/ {if ($2 == ""2"") print $7, $8, $9, $5 > ""msh_t.txt""}' ";

  command = [com_p,geometry,".msh ; "];
  command = [command,com_e,geometry,".msh ; "];
  command = [command,com_t,geometry,".msh"];
  
  system(command);

  ## 2d. Create PDE-tool like structure
  printf("Creating PDE-tool like mesh...\n");
  p   = load("msh_p.txt")'; # Mesh-points
  tmp = load("msh_e.txt")'; # Mesh surface-edges
  be  = zeros(7,columns(tmp));
  be([1,2,5],:) = tmp;
  t   = load("msh_t.txt")'; # Mesh tetrahedra

  ## 3. Remove hanging nodes
  printf("Check for hanging nodes...\n");
  nnodes = columns(p);
  in_msh = intersect( 1:nnodes , t(1:3,:) );
  if length(in_msh) != nnodes
    new_num(in_msh) = [1:length(in_msh)];
    t(1:3,:)        = new_num(t(1:3,:));
    be(1:2,:)       = new_num(be(1:2,:));
    p               = p(:,in_msh);
  endif

  ## 4. Set region numbers in edge structure
  printf("Setting region number in edge structure...\n");
  msh         = struct("p",p,"t",t,"e",be);
  msh.be      = MSH2Mtopprop(msh, "boundary");
  msh.e(6,:)  = msh.t(4,msh.be(1,:));
  jj          = find (sum(msh.be>0)==4);
  msh.e(7,jj) = msh.t(4,msh.be(3,jj)); 

  mesh = struct("p",p,"e",be,"t",t);

  ## 5. Delete temporary files
  printf("Deleting temporary files...\n");
  system(["rm -f msh_p.txt msh_e.txt msh_t.txt msh_s.txt *.msh"]);

endfunction

%!test
%! fid = fopen("circle.geo","w");
%! fprintf(fid,"Point(1) = {0, 0, 0, 1};\n");
%! fprintf(fid,"Point(2) = {1, 0, 0, 1};\n");
%! fprintf(fid,"Point(3) = {-1, 0, 0, 1};\n");
%! fprintf(fid,"Circle(1) = {3, 1, 2};\n");
%! fprintf(fid,"Circle(2) = {2, 1, 3};\n");
%! fprintf(fid,"Line Loop(4) = {2, 1};\n");
%! fprintf(fid,"Plane Surface(4) = {4};");
%! fclose(fid);
%! mesh = MSH2Mgmsh("circle","v",0);
%! system("rm circle.geo");
%! xidx = find(mesh.p(1,:) == 0);
%! yidx = find(mesh.p(2,:) == 0);
%! tmp  = intersect(xidx,yidx);
%! assert(isempty(tmp));