/* Copyright (C) 2013 Marco Vassallo
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <dolfin.h>
#include <octave/oct.h>
#include <octave/oct-map.h>

DEFUN_DLD (mshm_dolfin_write, args, ,"-*- texinfo -*-\n\
@deftypefn {Function File}\
mshm_dolfin_write (@var{mesh}, @var{mesh_name})\n\
Write a mesh to a dolfin .xml file.\n\
@itemize @bullet\n\
@item @var{mesh} is a PDE-tool like structure\n\
with matrix fields (p,e,t).\n\
@item The string @var{mesh_name} is an optional value specifying the output name.\n\
@end itemize\n\
@seealso{msh3m_structured_mesh, msh2m_structured_mesh, mshm_dolfin_read}\n\
@end deftypefn")
{
  int nargin = args.length ();
  octave_value_list retval;

  if (nargin < 1 || nargin > 2)
    print_usage ();
  else
    {
      Octave_map a (args(0).map_value ());
      if (! error_state)
        {
          std::string output_mesh;
          (nargin == 2) ? (output_mesh = args(1).string_value ()) : (output_mesh = "mesh");
          if (! error_state)
            {
              Array<octave_idx_type> p (a.contents (a.seek ("p"))(0).matrix_value ());
              Array<octave_idx_type> t (a.contents (a.seek ("t"))(0).matrix_value ());
              if (! error_state)
                {
                  dolfin::Mesh mesh;
                  std::size_t D = p.rows ();

                  if (D < 2 || D > 3)
                    error ("mshm_dolfin_write: only 2D or 3D meshes are supported");
                  else
                    {
                      dolfin::MeshEditor editor;
                      editor.open (mesh, D, D);
                      editor.init_vertices (p.cols ());
                      editor.init_cells (t.cols ());

                      if (D == 2)
                        {
                          for (uint i = 0; i < p.cols (); ++i)
                            editor.add_vertex (i, p.xelem (0, i), p.xelem (1, i));
      
                          for (uint i = 0; i < t.cols (); ++i)
                            editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1, t.xelem (2, i) - 1);
                        }

                      if (D == 3)
                        {
                          for (uint i = 0; i < p.cols (); ++i)
                            editor.add_vertex (i, p.xelem (0, i), p.xelem (1, i), p.xelem (2, i));
      
                          for (uint i = 0; i < t.cols (); ++i)
                            editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1, t.xelem (2, i) - 1, t.xelem (3, i) - 1);
                        }

                      editor.close ();

                      dolfin::File mesh_file (output_mesh + ".xml");
                      mesh_file << mesh;
                    }
                }
            }
        }
    }

  return retval;
}


/*
%!test
%! x = y = z = linspace (0, 1, 2);
%! msh = msh3m_structured_mesh (x, y, z, 1, [1 : 6]);
%! mshm_dolfin_write (msh, "msh");
%! msh = mshm_dolfin_read ("msh.xml");
%! p = [ 0   0   1   1   0   0   1   1
%!       0   1   0   1   0   1   0   1
%!       0   0   0   0   1   1   1   1];
%! assert (msh.p, p)
%! t = [   1   3   1   2   3   3
%!         2   5   3   3   6   4
%!         3   6   5   4   7   6
%!         6   7   6   6   8   8
%!         1   1   1   1   1   1];
%! assert (msh.t, t)
%! e = [1   1   5   3   1   1   2   2   6   3   4   3
%!      2   2   6   5   5   3   4   3   7   7   6   4
%!      6   3   7   7   6   5   6   4   8   8   8   8
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0
%!      0   0   0   0   0   0   0   0   0   0   0   0];
%! assert (msh.e, e)
*/
