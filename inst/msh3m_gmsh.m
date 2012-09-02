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
## msh3m_gmsh(@var{geometry}, @var{option}, @var{value}, ...) 
##
## Construct an unstructured tetrahedral 3D mesh making use of the free
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
## @seealso{msh3m_structured_mesh, msh2m_gmsh, msh2m_mesh_along_spline}
## @end deftypefn

function mesh = msh3m_gmsh (geometry, varargin)

  ## Check input
  ## Number of input
  if !(mod (nargin,2))
    warning("WRONG NUMBER OF INPUT.");
    print_usage;
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

  ## Generate mesh using Gmsh
  if (verbose)
    printf("\n");
    printf("Generating mesh...\n");
  endif
  system(["gmsh -format msh -3 -o " geometry ".msh" optstring " " geometry ".geo"]);
  
  if (verbose)
    printf("Processing gmsh data...\n");
  endif
  ## Points
  com_p   = "awk '/\\$Nodes/,/\\$EndNodes/ {print $2, $3, $4 > ""msh_p.txt""}' ";
  ## Surface edges
  com_e   = "awk '/\\$Elements/,/\\$EndElements/ {n=3+$3; if ($2 == ""2"") print $(n+1), $(n+2), $(n+3), $5 > ""msh_e.txt""}' ";
  ## Tetrahedra
  com_t   = "awk '/\\$Elements/,/\\$EndElements/ {n=3+$3; if ($2 == ""4"") print $(n+1), $(n+2), $(n+3), $(n+4), $5 > ""msh_t.txt""}' ";
  ## Side edges
  com_s   = "awk '/\\$Elements/,/\\$EndElements/ {n=3+$3; if ($2 == ""1"") print $(n+2), $(n+2), $5 > ""msh_s.txt""}' ";

  command = [com_p,geometry,".msh ; "];
  command = [command,com_e,geometry,".msh ; "];
  command = [command,com_t,geometry,".msh ; "];
  command = [command,com_s,geometry,".msh"];
  
  system(command);

  ## Create PDE-tool like structure
  if (verbose)
    printf("Creating PDE-tool like mesh...\n");
  endif
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
  if (verbose)
    printf("Check for hanging nodes...\n");
  endif
  nnodes = columns(p);
  in_msh = intersect( 1:nnodes , t(1:4,:) );
  if length(in_msh) != nnodes
    new_num(in_msh) = [1:length(in_msh)];
    t(1:4,:)        = new_num(t(1:4,:));
    be(1:3,:)       = new_num(be(1:3,:));
    p               = p(:,in_msh);
  endif

  mesh = struct("p",p,"s",s,"e",be,"t",t);
  
  if (verbose)
    printf("Deleting temporary files...\n");
  endif
  system(["rm -f msh_p.txt msh_e.txt msh_t.txt msh_s.txt *.msh"]);

endfunction