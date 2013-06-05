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

DEFUN_DLD(oct2xml,args,nargout,"oct2xml: oct2xml(my_mesh,output_name) \n convert a mesh from the [p,e,t] to the xml format and save it\n")
{
  int nargin = args.length();
  octave_value_list retval;

  if (nargin < 1 || nargin > 2)
    {
    print_usage ();
    }
  else
    {
    Octave_map a (args(0).map_value());
    if (! error_state)
      {
      std::string output_mesh;
      (nargin == 2)?(output_mesh = args(1).string_value()):(output_mesh = "mesh");

      Matrix p (a.contents (a.seek ("p"))(0).matrix_value());
      Matrix t (a.contents (a.seek ("t"))(0).matrix_value());
      if (! error_state)
	{
	std::size_t D = p.rows();
	std::size_t num_v = p.cols();
	std::size_t num_c = t.cols();
 
	dolfin::Mesh mesh;
	dolfin::MeshEditor editor;
 
	editor.open(mesh, D, D);
	editor.init_vertices(num_v);
	editor.init_cells(num_c);

	if( D == 2)
	  {
	  for(uint i = 0; i < num_v; ++i)
	    editor.add_vertex(i, p(0,i), p(1,i));
      
	  for(uint i = 0; i < num_c; ++i)
	    editor.add_cell(i, t(0,i)-1, t(1,i)-1, t(2,i)-1);
	  }

	if( D == 3)
	  {
	  for(uint i = 0; i < num_v;++i)
	    editor.add_vertex(i, p(0,i), p(1,i),p(2,i));
      
	  for(uint i = 0; i < num_c;++i)
	    editor.add_cell(i, t(0,i)-1 , t(1,i)-1, t(2,i)-1, t(3,i)-1 );
	  }

	editor.close();

	dolfin::File mesh_file(output_mesh+".xml");
	mesh_file << mesh;

	}
      }
    }
  return retval;
}


/*

*/
