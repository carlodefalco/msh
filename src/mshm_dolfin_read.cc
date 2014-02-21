/* Copyright (C) 2013 Marco Vassallo

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

DEFUN_DLD (mshm_dolfin_read, args, ,"-*- texinfo -*-\n\
@deftypefn {Function File} {[@var{mesh}]} = \
mshm_dolfin_read (@var{mesh_to_read}) \n\
Read a mesh from a dolfin .xml.gz file.\n\
The string @var{mesh_to_read} should be the name of the \
mesh file to be read.\n\
The output @var{mesh} is a PDE-tool like structure\n\
with matrix fields (p,e,t).\n\
@seealso{msh3m_structured_mesh, msh2m_structured_mesh, mshm_dolfin_write}\n\
@end deftypefn")
{
  octave_value_list retval;
#ifndef HAVE_DOLFIN_H
  error("mshm_dolfin_read: the msh package was built without support for dolfin (dolfin.h required)");
#else
  int nargin = args.length ();
  dim_vector dims;
  dims.resize (2);

  if (nargin != 1)
    print_usage ();
  else
    {
      std::string mesh_to_read = args(0).string_value ();
      if (! error_state)
        {
          boost::shared_ptr<dolfin::Mesh> mesh (new dolfin::Mesh (mesh_to_read));
          uint D = mesh->topology ().dim ();
          if (D < 2 || D > 3)
            error ("mshm_dolfin_read: only 2D or 3D meshes are supported");
          else
            {
              // matrix p
              std::size_t num_v = mesh->num_vertices ();
              Matrix p (D, num_v);
              std::copy (mesh->coordinates ().begin (),
                         mesh->coordinates ().end (),
                         p.fortran_vec ());

              // e has 7 rows in 2d, 10 rows in 3d
              mesh->init (D - 1, D);
              std::size_t num_f = mesh->num_facets ();
              dims(0) = D == 2 ? 7 : 10;
              dims(1) = num_f;
              Array<octave_idx_type> e (dims, 0);
              octave_idx_type *evec = e.fortran_vec ();
              uint D2 = D * D;
              octave_idx_type l = 0, m = 0;

              dolfin::MeshFunction <std::size_t> facet_domains (mesh, D - 1);
              bool empty = true;
              if (! mesh->domains ().is_empty ())
                  if (mesh->domains ().num_marked (D-1) != 0)
                    {
                      empty = false;
                      dolfin::MeshFunction<std::size_t> 
                        facet_domains_tmp (mesh, D - 1, mesh->domains ());
                      facet_domains = facet_domains_tmp;
                    }

              for (dolfin::FacetIterator f (*mesh); ! f.end (); ++f)
                {
                  if ((*f).exterior () == true)
                    {
                      l = 0;
                      for (dolfin::VertexIterator v (*f); ! v.end (); ++v, ++l)
                        e.xelem (l, m) = (*v).index () + 1;

                      if (! empty)
                        e.xelem (D2, m) = facet_domains[*f];

                      ++m;
                    }
                }

              dims(1) = m;
              e.resize (dims);

              for (octave_idx_type j = e.rows () - 2;
                   j < e.numel () - 2; j += e.rows ())
                evec[j] = 1;

              // t matrix
              dims(0) = D + 2;
              dims(1) = mesh->num_cells ();
              Array<octave_idx_type> t (dims, 1);
              std::vector<unsigned int> my_cells = mesh->cells ();
              std::size_t n = 0;

              empty = true;
              boost::shared_ptr<dolfin::Mesh> msh;
              dolfin::MeshFunction<std::size_t> cell_domains;
                if (! mesh->domains ().is_empty ())
                  if (mesh->domains ().num_marked (D) != 0)
                    {
                      empty = false;
                      dolfin::MeshFunction<std::size_t> 
                        cell_domains_tmp (mesh, D, mesh->domains ());
                      cell_domains = cell_domains_tmp;
                    }

              for (octave_idx_type j = 0; j < t.cols (); ++j)
                {
                  for (octave_idx_type i = 0; i < D + 1; ++i, ++n)
                    t.xelem (i, j) += my_cells[n];

                   if (! empty)
                     t.xelem (D + 1, j) = cell_domains[j];
                }

              octave_scalar_map a;
              a.setfield ("p", p);
              a.setfield ("e", e);
              a.setfield ("t", t);
              retval = octave_value (a);
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
