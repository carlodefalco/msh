## Copyright (C) 2013  Carlo de Falco
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

## -*- texinfo -*-
## @deftypefn {Function File} {} = msh3m_gmsh_write (@var{filename}, @var{msh})
## @seealso{msh2m_gmsh_write}
## @end deftypefn

function msh3m_gmsh_write (filename, msh, node_data, cell_data)
  
  if (! ((fid = fopen (filename, "w")) >= 0));
    error ("msh3m_gmsh_write: unable to open file %s for writing", filename);
  else
    ## file format string
    fprintf (fid, "$MeshFormat\n2.0 0 8\n$EndMeshFormat\n");

    ## node coordinates
    nnodes = columns (msh.p);
    fprintf (fid, "$Nodes\n%d\n", nnodes);
    p = [1:nnodes; msh.p];
    fprintf (fid, "%d %17.17g %17.17g %17.17g\n", p);
    fprintf (fid, "$EndNodes\n");

    ## elements
    number_of_tets = columns (msh.t);
    number_of_tri  = columns (msh.e);
    fprintf (fid, "$Elements\n%d\n", number_of_tets + number_of_tri);

    ## 3-node triangles 
    e = [1:number_of_tri;          ## element number
         2*ones(1, number_of_tri); ## element type, 2 = triangle
         3*ones(1, number_of_tri); ## number of tags
         zeros(1, number_of_tri);  ## first tag, physical entity: 0 = unspecified
         msh.e(10, :);             ## second tag, geometrical entity
         zeros(1, number_of_tri);  ## third tag, partition: 0 = unspecified
         msh.e(1:3, :)];           ## node number list
    fprintf (fid, "%d %d %d %d %d %d %d %d %d\n", e);

    ## 4-node tetrahedra
    t = [[(number_of_tri+1):(number_of_tets+number_of_tri)]; ## element number
         3*ones(1, number_of_tets);        ## element type, 3 = tetrahedron
         3*ones(1, number_of_tets);        ## number of tags
         zeros(1, number_of_tets);         ## first tag, physical entity: 0 = unspecified
         msh.t(5, :);                      ## first tag, geometrical entity
         zeros(1, number_of_tets);         ## third tag, partition: 0 = unspecified
         msh.t(1:4, :)];                   ## node number list
    fprintf (fid, "%d %d %d %d %d %d %d %d %d %d\n", t);
    fprintf(fid, "$EndElements\n");

    ## node data
    if (! isempty (node_data))
      for ii = 1:rows (node_data)
        fprintf (fid, "$NodeData\n")
        fprintf (fid, "%d\n", 1)                     ## number of string tags
        fprintf (fid, """%s""\n", node_data{ii, 1})  ## name of view
        fprintf (fid, "%d\n", 1)                     ## number of real tags
        fprintf (fid, "%g\n", 0.0)                   ## time
        fprintf (fid, "%d\n", 4)                     ## number of int tags
        fprintf (fid, "%d\n", [1, 1, nnodes, 0])     
        v = [1:nnodes; node_data{ii, 2}(:)'];
        fprintf (fid, "%d %g\n", v);
        fprintf (fid, "$EndNodeData\n");
      endfor
    endif
    fclose (fid);
  endif

endfunction
