## Copyright (C) 2006,2007,2008,2009,2010  Carlo de Falco, Massimiliano Culpo
##
## This file is part of:
##     MSH - Meshing Software Package for Octave
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
##  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>
##  author: Massimiliano Culpo <culpo _AT_ users.sourceforge.net>

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{mesh}]} = @
## msh2m_gmsh(@var{geometry},@var{option},@var{value},...)
##
## Construct an unstructured triangular 2D mesh making use of the free
## software gmsh.
##
## The compulsory argument @var{geometry} is the basename of the
## @code{*.geo} file to be meshed. 
##
## The optional arguments @var{option} and @var{value} identify
## respectively a gmsh option and its value. For more information
## regarding the possible option to pass, refer to gmsh manual or gmsh
## site @url{http://www.geuz.org/gmsh/}. 
##
## The returned value @var{mesh} is a PDE-tool like mesh structure.
##
## @seealso{msh2m_structured_mesh, msh3m_gmsh, msh2m_mesh_along_spline}
## @end deftypefn

function [mesh] = msh2m_gmsh(geometry,varargin)

  ## Check input
  if !mod(nargin,2) # Number of input parameters
    error("msh2m_gmsh: wrong number of input parameters.");
  endif
  ## FIXME: add input type check?

  ## Build mesh
  noptions  = (nargin - 1) / 2; # Number of passed options
  
  ## Construct system command string
  verbose   = 1;
  optstring = "";
  for ii = 1:noptions
    option    = varargin{2*(ii)-1};
    value     = varargin{2*ii};
    ## Check for verbose option
    if strcmp(option,"v")
      verbose = value;
    endif
    if !ischar(value)
      value = num2str(value);
    endif
    optstring = [optstring," -",option," ",value];
  endfor

  ## Invoke gmsh
  if (verbose)
    printf("\n");
    printf("Generating mesh...\n");
  endif
  system(["gmsh -format msh -2 -o " geometry ".msh" optstring " " geometry ".geo"]);

  ## Build structure fields
  if (verbose)
    printf("Processing gmsh data...\n");
  endif
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

  ## Create PDE-tool like structure
  if (verbose)
    printf("Creating PDE-tool like mesh...\n");
  endif
  p   = load("msh_p.txt")'; # Mesh-points
  tmp = load("msh_e.txt")'; # Mesh surface-edges
  be  = zeros(7,columns(tmp));
  be([1,2,5],:) = tmp;
  t   = load("msh_t.txt")'; # Mesh tetrahedra

  ## Remove hanging nodes
  if (verbose)
    printf("Check for hanging nodes...\n");
  endif
  nnodes = columns(p);
  in_msh = intersect( 1:nnodes , t(1:3,:) );
  if length(in_msh) != nnodes
    new_num(in_msh) = [1:length(in_msh)];
    t(1:3,:)        = new_num(t(1:3,:));
    be(1:2,:)       = new_num(be(1:2,:));
    p               = p(:,in_msh);
  endif

  ## Set region numbers in edge structure
  if (verbose)
    printf("Setting region number in edge structure...\n");
  endif
  msh         = struct("p",p,"t",t,"e",be);
  msh.be      = msh2m_topological_properties(msh, "boundary");
  msh.e(6,:)  = msh.t(4,msh.be(1,:));
  jj          = find (sum(msh.be>0)==4);
  msh.e(7,jj) = msh.t(4,msh.be(3,jj)); 

  mesh = struct("p",p,"e",be,"t",t);

  ## Delete temporary files
  if (verbose)
    printf("Deleting temporary files...\n");
  endif
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
%! mesh = msh2m_gmsh("circle","v",0);
%! system("rm circle.geo");
%! nnodest = length(unique(mesh.t));
%! nnodesp = columns(mesh.p);
%! assert(nnodest,nnodesp);