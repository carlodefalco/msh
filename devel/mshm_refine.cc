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

DEFUN_DLD (mshm_refine, args, ,
	   "-*- texinfo -*-\n\
@deftypefn {Function File} @var{refined_mesh} @\n\
mshm_refine(@var{mesh})\n\
\n\
Uniformly refine a mesh.\n\
\n\
@itemize @bullet\n\
@item @var{mesh} is the mesh you want to be refined \n\
@end itemize\n\
\n\
WARNING: the @var{refined_mesh} is a new mesh, which doesn't share any \n\
information with the original @var{mesh}.  \n\
For example, border labels could be different.   \n\
@seealso{msh3m_structured_mesh, msh2m_structured_mesh}\n\
@end deftypefn")
{
  int nargin = args.length ();
  octave_value_list retval;

  if (nargin != 1)
    {
      print_usage ();
    }
  else
    {
      Octave_map a (args(0).map_value ());
      if (! error_state)
	{
	  Matrix p (a.contents (a.seek ("p"))(0).matrix_value ());
	  Matrix t (a.contents (a.seek ("t"))(0).matrix_value ());
	  if (! error_state)
	    {
	      dolfin::Mesh mesh;
	      std::size_t D = p.rows ();
	      dolfin::MeshEditor editor;
	      editor.open (mesh, D, D);

	      std::size_t num_v = p.cols ();
	      editor.init_vertices (num_v);

	      std::size_t num_c = t.cols ();
	      editor.init_cells (num_c);

	      if (D == 2)
		{
		  for (uint i = 0; i < num_v; ++i)
		    editor.add_vertex (i, p(0,i), p(1,i));
      
		  for (uint i = 0; i < num_c; ++i)
		    editor.add_cell (i, t(0,i)-1, t(1,i)-1, t(2,i)-1);
		}

	      if (D == 3)
		{
		  for (uint i = 0; i < num_v; ++i)
		    editor.add_vertex (i, p(0,i), p(1,i), p(2,i));
      
		  for (uint i = 0; i < num_c; ++i)
		    editor.add_cell (i, t(0,i)-1, t(1,i)-1, t(2,i)-1, t(3,i)-1);
		}

	      editor.close ();

	      dolfin::Mesh refined_mesh;
	      refined_mesh = dolfin::refine(mesh);

	      std::size_t num_v_r = refined_mesh.num_vertices ();
	      Matrix p_r (D, num_v_r);
	      std::vector<double> my_coord_r = refined_mesh.coordinates ();
	      std::size_t n = 0;

	      for (octave_idx_type j=0; j < num_v_r; ++j)
		for (octave_idx_type i =0; i < D; ++i,++n)
		  p_r(i,j) = my_coord_r[n];
	    

	      refined_mesh.init (D - 1, D);
	      std::size_t num_f_r = refined_mesh.num_facets ();
	      Matrix e_r;
	      if(D == 2)
		e_r.resize (7, num_f_r);
	      if(D == 3)
		e_r.resize (10, num_f_r);

	      octave_idx_type l=0;
	      octave_idx_type m=0;

	      for (dolfin::FacetIterator f (refined_mesh); !f.end (); ++f)
		{
		  if ( (*f).exterior () == true)
		    {
		      l = 0;
		      for (dolfin::VertexIterator v (*f); !v.end (); ++v,++l)
			{ 
			  e_r(l,m) = (*v).index () + 1 ;
			}
		      ++m;
		    }   
		}

	      if (D == 2) 
		{
		  e_r.resize (7, m);
		  e_r.fill (1, 6, 0, 6, m-1);
		}

	      if (D == 3) 
		{
		  e_r.resize (10, m);
		  e_r.fill (1, 9, 0, 9, m-1);
		}
  
	  
	      std::size_t num_c_r = refined_mesh.num_cells ();
	      Matrix t_r (D+2, num_c_r);
	      t_r.fill (1.0);
	      std::vector<unsigned int> my_cells_r = refined_mesh.cells ();
	      n = 0;

	      for (octave_idx_type j=0; j < num_c_r; ++j)
		for (octave_idx_type i =0; i < D+1; ++i,++n)
		  t_r(i,j) += my_cells_r[n]; 

	      Octave_map a_r;
	      a_r.assign ("p", octave_value (p_r));
	      a_r.assign ("e", octave_value (e_r));
	      a_r.assign ("t", octave_value (t_r));
	      retval = octave_value (a_r);


	    }
	}	
    }
  return retval;
}


/*
%!demo
%! msh = msh2m_structured_mesh(linspace(0,1,3),linspace(0,1,3),1,[1:4]);
%! refined_mesh = mshm_refine(msh);
%! clf;
%! subplot(1,2,1);
%! msh2p_mesh(msh);
%! subplot(1,2,2);
%! msh2p_mesh(refined_mesh);
*/
