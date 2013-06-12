
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

DEFUN_DLD (mshm_refine, args, ,"-*- texinfo -*-\n\
@deftypefn {Function File} {[@var{refined_mesh}]} = \
mshm_refine (@var{mesh},@var{cell_marker}) \n\
Refine a mesh\n\
@itemize @bullet \n\
@item The @var{mesh} is a PDE-tool like structures with matrix field (p,e,t).\n\
@item The optional argument @var{cell_marker} is a list\n\
containing the number of the cells you want to refine.\n\
By default a uniform refinement is applied.\n\
@end itemize\n\
The output @var{refined_mesh} is a refined mesh with\n\
the same structure as @var{mesh}\n\
@seealso{msh3m_structured_mesh, msh2m_structured_mesh}\n\
@end deftypefn")
{
  int nargin = args.length ();
  octave_value_list retval;
  dim_vector dims;
  dims.resize (2);

  if (nargin < 1 || nargin > 2)
    print_usage ();
  else
    {
      Octave_map a (args(0).map_value ());
      if (! error_state)
        {
          Array<double> p (a.contents (a.seek ("p"))(0).matrix_value ());
          Array<octave_idx_type> t (a.contents (a.seek ("t"))(0).matrix_value ());
          if (! error_state)
            {
              dolfin::Mesh mesh;
              std::size_t D = p.rows ();

              if (D < 2 || D > 3)
                error ("mshm_refine: only 2D or 3D meshes are supported");
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
                        editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1,
                                            t.xelem (2, i) - 1);
                    }

                  if (D == 3)
                    {
                      for (uint i = 0; i < p.cols (); ++i)
                        editor.add_vertex (i, p.xelem (0, i), p.xelem (1, i), p.xelem (2, i));

                      for (uint i = 0; i < t.cols (); ++i)
                        editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1,
                                            t.xelem (2, i) - 1, t.xelem (3, i) - 1);
                    }

                  editor.close ();

                  if (nargin == 2)
                    { 
                      Array<octave_idx_type> cell_idx (args(1).array_value ());
                      if (! error_state)
                        {
                          dolfin::CellFunction<bool> cell_markers (mesh);
                          cell_markers.set_all (false);

                          for (octave_idx_type i = 0; i < cell_idx.length (); ++i)
                            cell_markers.set_value (cell_idx.xelem (i) - 1, true);

                          mesh = dolfin::refine (mesh, cell_markers);
                        }
                    }
                  else
                    mesh = dolfin::refine(mesh);


                  std::size_t num_v = mesh.num_vertices ();
                  dims(0) = D;
                  dims(1) = num_v;
                  p.resize (dims);
                  std::copy (mesh.coordinates ().begin (),
                             mesh.coordinates ().end (),
                             p.fortran_vec ());
                
                  mesh.init (D - 1, D);
                  std::size_t num_f = mesh.num_facets ();
                
                  // e has 7 rows in 2d, 10 rows in 3d
                  dims(0) = D == 2 ? 7 : 10; 
                  dims(1) = num_f;
                  Array<octave_idx_type> e (dims, 0);
                  octave_idx_type *evec = e.fortran_vec ();

                  octave_idx_type l = 0, m = 0;

                  for (dolfin::FacetIterator f (mesh); ! f.end (); ++f)
                    {
                      if ((*f).exterior () == true)
                        {
                          l = 0;
                          for (dolfin::VertexIterator v (*f); ! v.end (); ++v, ++l) 
                            e.xelem (l, m) = (*v).index () + 1 ;
                          ++m;
                        }
                    }

                  dims(1) = m;
                  e.resize (dims);

                  for (octave_idx_type j = e.rows () - 2; 
                       j < e.numel () - 2; j += e.rows ())
                    evec[j] = 1;

                  dims(0) = D + 2;
                  dims(1) =  mesh.num_cells ();
                  t.clear ();
                  t.resize (dims, 1);

                  std::vector<unsigned int> my_cells = mesh.cells ();
                  std::size_t n = 0;

                  for (octave_idx_type j = 0; j < t.cols (); ++j)
                    for (octave_idx_type i = 0; i < D + 1; ++i, ++n)
                      t.xelem (i, j) += my_cells[n]; 

                  Octave_map a;
                  a.assign ("p", p);
                  a.assign ("e", e);
                  a.assign ("t", t);
                  retval = octave_value (a);


                }
            }
        }
    }

  return retval;
}


/*
%!demo
%! msh = msh2m_structured_mesh (linspace (0, 1, 4), linspace (0, 1, 4), 1, [1 : 4]);
%! partially_refined_mesh = mshm_refine (msh,[1:3]);
%! uniformly_refined_mesh = mshm_refine (msh);
%! clf;
%! subplot (1, 3, 1);
%! msh2p_mesh (msh);
%! title ('original mesh');
%! subplot (1, 3, 2);
%! msh2p_mesh (partially_refined_mesh);
%! title ('partially refined mesh');
%! subplot (1, 3, 3);
%! msh2p_mesh (uniformly_refined_mesh);
%! title ('uniformly refined mesh');
*/
