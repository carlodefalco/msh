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

#ifdef HAVE_DOLFIN_H
#include <dolfin.h>
#endif
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
#ifndef HAVE_DOLFIN_H
  error("mshm_dolfn_write: the msh package was built without support for dolfin (dolfin.h required)");
#else

  if (nargin < 1 || nargin > 2)
    print_usage ();
  else
    {
      octave_scalar_map a = args(0).scalar_map_value ();
      std::string output_mesh;

      output_mesh = "mesh";
      if (nargin == 2)
        output_mesh = args(1).string_value ();

      Array<double> p = a.contents ("p").matrix_value ();
      Array<octave_idx_type> t = a.contents ("t").matrix_value ();
      Array<octave_idx_type> e = a.contents ("e").matrix_value ();
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
                    editor.add_cell (i, t.xelem (0, i) - 1,
                                     t.xelem (1, i) - 1, t.xelem (2, i) - 1);
                }

              if (D == 3)
                {
                  for (uint i = 0; i < p.cols (); ++i)
                    editor.add_vertex (i, p.xelem (0, i),
                                       p.xelem (1, i), p.xelem (2, i));

                  for (uint i = 0; i < t.cols (); ++i)
                    editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1,
                                     t.xelem (2, i) - 1, t.xelem (3, i) - 1);
                }

              editor.close ();

              // store information associated with e
              mesh.init (D - 1);
              dolfin::MeshValueCollection<uint> facet (D - 1);
              std::size_t num_side_edges = e.cols ();

              if (D == 2)
                {
                  for (uint i = 0; i < num_side_edges; ++i)
                    {
                      dolfin::Vertex v (mesh, e.xelem (0, i) - 1);
                      for (dolfin::FacetIterator f (v); ! f.end (); ++f)
                        {
                          if ((*f).entities(0)[0] == e.xelem (0, i) - 1
                              && (*f).entities(0)[1] == e.xelem (1, i) - 1
                              || (*f).entities(0)[0] == e.xelem (1, i) - 1
                              && (*f).entities(0)[1] == e.xelem (0, i) - 1)
                            {
                              facet.set_value ((*f).index (), e.xelem (4, i), mesh);
                              break;
                            }
                        }
                    }
                }



              if (D == 3)
                {
                  for (uint i = 0; i < num_side_edges; ++i)
                    {
                      dolfin::Vertex v (mesh, e.xelem (0, i) - 1);
                      for (dolfin::FacetIterator f (v); ! f.end (); ++f)
                        {
                          if ((*f).entities(0)[0] == e(0, i) - 1
                              && (*f).entities(0)[1] == e.xelem (1, i) - 1
                              && (*f).entities(0)[2] == e.xelem (2, i) - 1
                              || (*f).entities(0)[0] == e.xelem (0, i) - 1
                              && (*f).entities(0)[1] == e.xelem (2, i) - 1
                              && (*f).entities(0)[2] == e.xelem (1, i) - 1
                              || (*f).entities(0)[0] == e.xelem (1, i) - 1
                              && (*f).entities(0)[1] == e.xelem (0, i) - 1
                              && (*f).entities(0)[2] == e.xelem (2, i) - 1
                              || (*f).entities(0)[0] == e.xelem (1, i) - 1
                              && (*f).entities(0)[1] == e.xelem (2, i) - 1
                              && (*f).entities(0)[2] == e.xelem (0, i) - 1
                              || (*f).entities(0)[0] == e.xelem (2, i) - 1
                              && (*f).entities(0)[1] == e.xelem (0, i) - 1
                              && (*f).entities(0)[2] == e.xelem (1, i) - 1
                              || (*f).entities(0)[0] == e.xelem (2, i) - 1
                              && (*f).entities(0)[1] == e.xelem (1, i) - 1
                              && (*f).entities(0)[2] == e.xelem (0, i) - 1)
                            {
                              facet.set_value ((*f).index (), e.xelem (9, i), mesh);
                              break;
                            }
                        }
                    }
                }

              *(mesh.domains ().markers (D - 1)) = facet;

              // store information associated with t
              dolfin::MeshValueCollection<uint> cell (D);
              std::size_t num_cells = t.cols ();

              if (D == 2)
                {
                  for (uint i = 0; i < num_cells; ++i)
                    {
                      dolfin::Vertex v (mesh, t.xelem (0, i) - 1);
                      for (dolfin::CellIterator f (v); ! f.end (); ++f)
                        {
                          if ((*f).entities(0)[0] == t.xelem (0, i) - 1
                              && (*f).entities(0)[1] == t.xelem (1, i) - 1
                              && (*f).entities(0)[2] == t.xelem (2, i) - 1
                              || (*f).entities(0)[0] == t.xelem (0, i) - 1
                              && (*f).entities(0)[1] == t.xelem (2, i) - 1
                              && (*f).entities(0)[2] == t.xelem (1, i) - 1
                              || (*f).entities(0)[0] == t.xelem (1, i) - 1
                              && (*f).entities(0)[1] == t.xelem (0, i) - 1
                              && (*f).entities(0)[2] == t.xelem (2, i) - 1
                              || (*f).entities(0)[0] == t.xelem (1, i) - 1
                              && (*f).entities(0)[1] == t.xelem (2, i) - 1
                              && (*f).entities(0)[2] == t.xelem (0, i) - 1
                              || (*f).entities(0)[0] == t.xelem (2, i) - 1
                              && (*f).entities(0)[1] == t.xelem (0, i) - 1
                              && (*f).entities(0)[2] == t.xelem (1, i) - 1
                              || (*f).entities(0)[0] == t.xelem (2, i) - 1
                              && (*f).entities(0)[1] == t.xelem (1, i) - 1
                              && (*f).entities(0)[2] == t.xelem (0, i) - 1)
                            {
                              cell.set_value ((*f).index (), t.xelem (3, i), mesh);
                              break;
                            }
                        }
                    }
                }



              if (D == 3)
                {
                  for (uint i = 0; i < num_cells; ++i)
                    {
                      dolfin::Vertex v (mesh, t.xelem (0, i) - 1);
                      for (dolfin::CellIterator f (v); ! f.end (); ++f)
                        {
                          if ((*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1
                               || (*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (0, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1

                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1
                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (1, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1

                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1
                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (3, i) - 1
                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (3, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1
                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1
                               || (*f).entities(0)[0] == t.xelem (2, i) - 1
                               && (*f).entities(0)[1] == t.xelem (3, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1

                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (0, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1
                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (2, i) - 1
                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (1, i) - 1
                               && (*f).entities(0)[2] == t.xelem (2, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1
                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (0, i) - 1
                               && (*f).entities(0)[3] == t.xelem (1, i) - 1
                               || (*f).entities(0)[0] == t.xelem (3, i) - 1
                               && (*f).entities(0)[1] == t.xelem (2, i) - 1
                               && (*f).entities(0)[2] == t.xelem (1, i) - 1
                               && (*f).entities(0)[3] == t.xelem (0, i) - 1)
                            {
                              cell.set_value ((*f).index (), t.xelem (4, i), mesh);
                              break;
                            }
                        }
                    }
                }

              *(mesh.domains ().markers (D)) = cell;

              dolfin::File mesh_file (output_mesh + ".xml");
              mesh_file << mesh;
            }
        }
    }
#endif
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
%!      1   5   6   3   1   3   4   5   6   2   4   2];
%! assert (msh.e, e)
*/
