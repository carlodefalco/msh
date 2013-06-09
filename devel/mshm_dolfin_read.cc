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

DEFUN_DLD (mshm_dolfin_read, args, ,
"-*- texinfo -*-\n\
@deftypefn {Function File} {[@var{mesh}]} = @\n\
mshm_dolfin_read(@var{mesh_to_read})\n\
\n\
Read a dolfin mesh in .xml format and import it to Octave.\n\
\n\
@itemize @bullet\n\
@item @var{mesh_to_read} is the name of the mesh you want to read.\n\
@end itemize\n\
\n\
The returned value @var{mesh} is a PDE-tool like structure\n\
composed of the matrices (p,e,t).\n\
\n\
@seealso{msh3m_structured_mesh, msh2m_structured_mesh, mshm_dolfin_write}\n\
@end deftypefn")
{
  int nargin = args.length ();
  octave_value_list retval;

  if (nargin != 1)
    print_usage ();
  else
    {
      std::string mesh_to_read = args(0).string_value ();
      if (! error_state)
	{
	  dolfin::Mesh mesh (mesh_to_read);
	  uint D = mesh.topology ().dim ();

	  std::size_t num_v = mesh.num_vertices ();
	  Matrix p (D, num_v);
	  std::vector<double> my_coord = mesh.coordinates ();
	  std::size_t n = 0;

	  for (octave_idx_type j=0; j < num_v; ++j)
	    for (octave_idx_type i =0; i < D; ++i,++n)
	      p(i,j) = my_coord[n];
	    

	  mesh.init (D - 1, D);
	  std::size_t num_f = mesh.num_facets ();
	  Matrix e;
	  if(D == 2)
	    e.resize (7, num_f);
	  if(D == 3)
	    e.resize (10, num_f);

	  octave_idx_type l=0;
	  octave_idx_type m=0;

	  for (dolfin::FacetIterator f (mesh); !f.end (); ++f)
	    {
	      if ( (*f).exterior () == true)
		{
		  l = 0;
		  for (dolfin::VertexIterator v (*f); !v.end (); ++v,++l)
		    { 
		      e(l,m) = (*v).index () + 1 ;
		    }
		  ++m;
		}   
	    }

	  if (D == 2) 
	    {
	      e.resize (7, m);
	      e.fill (1, 6, 0, 6, m-1);
	    }

	  if (D == 3) 
	    {
	      e.resize (10, m);
	      e.fill (1, 9, 0, 9, m-1);
	    }
  
	  
	  std::size_t num_c = mesh.num_cells ();
	  Matrix t (D+2, num_c);
	  t.fill (1.0);
	  std::vector<unsigned int> my_cells = mesh.cells ();
	  n = 0;

	  for (octave_idx_type j=0; j < num_c; ++j)
	    for (octave_idx_type i =0; i < D+1; ++i,++n)
	      t(i,j) += my_cells[n]; 

	  Octave_map a;
	  a.assign ("p", octave_value (p));
	  a.assign ("e", octave_value (e));
	  a.assign ("t", octave_value (t));
	  retval = octave_value (a);
	}
    }
  return retval;
}


/*
%!demo
%! mesh = mshm_dolfin_read("dolfin_fine.xml.gz" );
%! msh2p_mesh(mesh);

*/
