/* Copyright (C) 2013-14 Marco Vassallo

   This file is part of:
   MSH - Meshing Software Package for Octave

   MSH is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   MSH is distributed in the hope that it will be useful,
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
#include <algorithm>

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
#ifndef HAVE_DOLFIN_H
  error("mshm_refine: the msh package was built without support for dolfin (dolfin.h required)");
#else
  dim_vector dims;
  dims.resize (2);

  if (nargin < 1 || nargin > 2)
    print_usage ();
  else
    {
      octave_scalar_map a = args(0).scalar_map_value ();
      Array<double> p = a.contents ("p").matrix_value ();
      Array<octave_idx_type> t = a.contents ("t").matrix_value ();
      Array<octave_idx_type> e = a.contents ("e").matrix_value ();
      Array<octave_idx_type> cell_idx;

      if (nargin == 2)
        cell_idx = args(1).array_value ();
      
      if (! error_state)
        {
           int min = *std::min_element (cell_idx.fortran_vec (),
                                        cell_idx.fortran_vec ()
                                        + cell_idx.length ());
           int max = *std::max_element (cell_idx.fortran_vec (),
                                        cell_idx.fortran_vec ()
                                        + cell_idx.length ());

           if (nargin == 2 && (min < 1 || max > t.cols ()))
             error ("mshm_refine: cell index out of bounds");
          else
            {
               boost::shared_ptr<dolfin::Mesh> mesh (new dolfin::Mesh ());
               std::size_t D = p.rows ();

               if (D < 2 || D > 3)
                 error ("mshm_refine: only 2D or 3D meshes are supported");
               else
                 {
                   dolfin::MeshEditor editor;
                   editor.open (*mesh, D, D);
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
                       for (unsigned int i = 0; i < p.cols (); ++i)
                         editor.add_vertex (i, p.xelem (0, i),
                                            p.xelem (1, i), p.xelem (2, i));

                       for (unsigned int i = 0; i < t.cols (); ++i)
                         editor.add_cell (i, t.xelem (0, i) - 1, t.xelem (1, i) - 1,
                                          t.xelem (2, i) - 1, t.xelem (3, i) - 1);
                     }

                   editor.close ();

                   // store information associated with e
                   mesh->init (D - 1);
                   std::size_t num_side_edges = e.cols ();

                   if (D == 2)
                     {
                       for (uint i = 0; i < num_side_edges; ++i)
                         {
                           dolfin::Vertex v (*mesh, e.xelem (0, i) - 1);
                           for (dolfin::FacetIterator f (v); ! f.end (); ++f)
                             {
                               if ((*f).entities(0)[0] == e.xelem (0, i) - 1
                                   && (*f).entities(0)[1] == e.xelem (1, i) - 1
                                   || (*f).entities(0)[0] == e.xelem (1, i) - 1
                                   && (*f).entities(0)[1] == e.xelem (0, i) - 1)
                                 {
                                    std::pair <std::size_t, std::size_t>
                                      idxvl ((*f).index (), e.xelem (4, i));
                                    mesh->domains ().set_marker (idxvl, D - 1);
                                    break;
                                 }
                             }
                         }
                     }

                   if (D == 3)
                     {
                       for (uint i = 0; i < num_side_edges; ++i)
                         {
                           dolfin::Vertex v (*mesh, e.xelem (0, i) - 1);
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
                                    std::pair <std::size_t, std::size_t> 
                                      idxvl ((*f).index (), e.xelem (9, i));
                                    mesh->domains ().set_marker (idxvl, D - 1);
                                    break;
                                 }
                             }
                         }
                     }

                   // store information associated with t
                   std::size_t num_cells = t.cols ();

                   if (D == 2)
                     {
                       for (uint i = 0; i < num_cells; ++i)
                         {
                           dolfin::Vertex v (*mesh, t.xelem (0, i) - 1);
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
                                    std::pair <std::size_t, std::size_t>
                                      idxvl ((*f).index (), t.xelem (3, i));
                                    mesh->domains ().set_marker (idxvl, D);
                                    break;
                                 }
                             }
                         }
                     }

                   if (D == 3)
                     {
                       for (uint i = 0; i < num_cells; ++i)
                         {
                           dolfin::Vertex v (*mesh, t.xelem (0, i) - 1);
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
                                    std::pair <std::size_t, std::size_t>
                                      idxvl ((*f).index (), t.xelem (4, i));
                                    mesh->domains ().set_marker (idxvl, D);
                                    break;
                                 }
                             }
                         }
                     }

                   dolfin::MeshFunction<std::size_t> 
                     cell (mesh, D, mesh->domains ());
                   dolfin::MeshFunction<std::size_t> 
                     facet (mesh, D - 1, mesh->domains ());

                   dolfin::CellFunction<bool> cell_markers (*mesh);
                   if (nargin == 2)
                     {
                       cell_markers.set_all (false);
                       for (octave_idx_type i = 0; i < cell_idx.length (); ++i)
                         cell_markers.set_value (cell_idx.xelem (i) - 1, true);
                     }
                   else
                     cell_markers.set_all (true);

                   boost::shared_ptr<dolfin::Mesh> r_mesh(new dolfin::Mesh ());
                   dolfin::refine (*r_mesh, *mesh, cell_markers);

                   std::size_t num_v = (*r_mesh).num_vertices ();
                   dims(0) = D;
                   dims(1) = num_v;
                   p.resize (dims);
                   std::copy ((*r_mesh).coordinates ().begin (),
                              (*r_mesh).coordinates ().end (),
                              p.fortran_vec ());

                   // e has 7 rows in 2d, 10 rows in 3d
                   (*r_mesh).init (D - 1, D);
                   std::size_t num_f = (*r_mesh).num_facets ();
                   dims(0) = D == 2 ? 7 : 10;
                   dims(1) = num_f;
                   e.clear ();
                   e.resize (dims, 0);
                   octave_idx_type *evec = e.fortran_vec ();
                   uint D2 = D * D;
                   octave_idx_type l = 0, m = 0;

                   dolfin::MeshFunction<std::size_t> r_facet (*r_mesh, D - 1);
                   r_facet = dolfin::adapt (facet, r_mesh);

                   for (dolfin::FacetIterator f (*r_mesh); ! f.end (); ++f)
                     {
                       if ((*f).exterior () == true)
                         {
                            l = 0;
                            for (dolfin::VertexIterator v (*f); ! v.end (); ++v, ++l)
                              e.xelem (l, m) = (*v).index () + 1 ;

                            e.xelem (D2, m) = r_facet[*f];
                            ++m;
                         }
                     }

                   dims(1) = m;
                   e.resize (dims);

                   for (octave_idx_type j = e.rows () - 2;
                        j < e.numel () - 2; j += e.rows ())
                     evec[j] = 1;

                   dims(0) = D + 2;
                   dims(1) = (*r_mesh).num_cells ();
                   t.clear ();
                   t.resize (dims, 1);
                   std::vector<unsigned int> my_cells = (*r_mesh).cells ();
                   std::size_t n = 0;

                   dolfin::MeshFunction<std::size_t> r_cell (*r_mesh, D);
                   r_cell = dolfin::adapt (cell, r_mesh);

                   for (octave_idx_type j = 0; j < t.cols (); ++j)
                     {
                       for (octave_idx_type i = 0; i < D + 1; ++i, ++n)
                         t.xelem (i, j) += my_cells[n];

                       t.xelem (D + 1, j) = r_cell[j];
                     }

                   a.setfield ("p", p);
                   a.setfield ("e", e);
                   a.setfield ("t", t);
                   retval = octave_value (a);
                }
            }
        }
    }
#endif
  return retval;
}

/*
%!demo
%! # Create a uniform mesh 
%! msh = msh2m_structured_mesh (linspace (0, 1, 4), linspace (0, 1, 4), 1, [1 : 4]);
%! # Refine it only on cells from 1 to 3
%! partially_refined_mesh = mshm_refine (msh,[1:3]);
%! # Refine the original mesh uniformly
%! uniformly_refined_mesh = mshm_refine (msh);
%!
%! # plot the result 
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
/*
%!test
%! x = y = linspace (0, 1, 2);
%! msh = msh2m_structured_mesh (x, y, 1, [1 : 4]);
%! msh.t (4, 2) = 2;
%! msh_r = mshm_refine (msh);
%! msh_rr = mshm_refine (msh_r);
%! p = [ 0.00000   0.00000   1.00000   1.00000   0.50000   0.50000   1.00000   0.00000   0.50000
%!       0.00000   1.00000   0.00000   1.00000   0.50000   0.00000   0.50000   0.50000   1.00000];
%! assert (msh_rr.p, p)
%! t = [ 1   3   3   4   1   2   2   4
%!       5   5   5   5   5   5   5   5
%!       6   6   7   7   8   8   9   9
%!       1   1   1   1   2   2   2   2];
%! assert (msh_rr.t, t)
%! e =[  1   3   3   4   1   2   2   4
%!       6   6   7   7   8   8   9   9
%!       0   0   0   0   0   0   0   0
%!       0   0   0   0   0   0   0   0
%!       1   1   2   2   4   4   3   3
%!       0   0   0   0   0   0   0   0
%!       0   0   0   0   0   0   0   0];
%! assert (msh_rr.e, e)
*/
