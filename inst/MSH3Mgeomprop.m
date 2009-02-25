## Copyright (C) 2007,2008  Carlo de Falco, Massimiliano Culpo
##
## This file is part of 
##
##                   MSH - Meshing Software Package for Octave
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
##
##  AUTHORS:
##  Carlo de Falco
##  Dublin City University
##  School of Mathemetical Sciences
##  Ireland
##
##  Culpo Massimiliano
##  Bergische Universitaet Wuppertal
##  Fachbereich C - Mathematik und Naturwissenschaften
##  Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
##  D-42119 Wuppertal, Germany

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{varargout}]} = MSH3Mgeomprop(@var{mesh},[@var{string1},@var{string2},...])
##
## Computes geometrical properties of the specified mesh
##
## Input:
## @itemize @minus
## @item @var{mesh}: standard PDEtool-like mesh, with field "p", "e", "t".
## @item @var{string1}, @var{string2},...: identifier of the property to compute. Possible choices are listed below.
## @itemize @bullet
## @item "wjacdet" : weighted determinant of jacobian trasformation to
## reference tetrahedron
## @item "shg": gradient of the P1 shape functions for BIM method
## @item "shp": value of the P1 shape functions for BIM method
## @end itemize
## @end itemize 
##
## The output will contain the geometrical properties requested in the input in the same order specified in the function call
## @end deftypefn


function [varargout] = MSH3Mgeomprop(imesh,varargin)

  ## check if varargin is composed of strings as expected
  if iscellstr(varargin) == 0
    warning("Unexpected input. See help for more information.");
    print_usage;
  endif
  
  ## extract tetrahedra node coordinates
  x1 = imesh.p(1,imesh.t(1,:));
  y1 = imesh.p(2,imesh.t(1,:));
  z1 = imesh.p(3,imesh.t(1,:));
  x2 = imesh.p(1,imesh.t(2,:));
  y2 = imesh.p(2,imesh.t(2,:));
  z2 = imesh.p(3,imesh.t(2,:));
  x3 = imesh.p(1,imesh.t(3,:));
  y3 = imesh.p(2,imesh.t(3,:));
  z3 = imesh.p(3,imesh.t(3,:));
  x4 = imesh.p(1,imesh.t(4,:));
  y4 = imesh.p(2,imesh.t(4,:));
  z4 = imesh.p(3,imesh.t(4,:));

  for nn = 1:length(varargin)
    
    request = varargin{nn};
    switch request
      case "wjacdet"
	b = wjacdet(x1,y1,z1,\
		    x2,y2,z2,\
		    x3,y3,z3,\
		    x4,y4,z4);
	varargout{nn} = b;
        clear b
      case "area"
	tmp = wjacdet(x1,y1,z1,\
		      x2,y2,z2,\
		      x3,y3,z3,\
		      x4,y4,z4);
	b   = sum(tmp,1);
	varargout{nn} = b;
      case "shg"
	b = shg(x1,y1,z1,\
		x2,y2,z2,\
		x3,y3,z3,\
		x4,y4,z4);
	varargout{nn} = b;
        clear b
      case "shp"
	varargout{nn} = eye(4);
      otherwise
	warning("Unexpected value in passed string. Empty vector passed as output.")
	varargout{nn} = [];
    endswitch
    
  endfor

endfunction

function [b] = wjacdet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
  ## COMPUTE WEIGHTED YACOBIAN DETERMINANT
  
  ## weight
  weight = [1/4 1/4 1/4 1/4]';

  Nb2 = y1.*(z3-z4) + y3.*(z4-z1) + y4.*(z1-z3);
  Nb3 = y1.*(z4-z2) + y2.*(z1-z4) + y4.*(z2-z1);
  Nb4 = y1.*(z2-z3) + y2.*(z3-z1) + y3.*(z1-z2);
  
  ## Determinant of the Jacobian of the 
  ## transformation from the base tetrahedron
  ## to the tetrahedron K  
  detJ = (x2-x1).*Nb2 +(x3-x1).*Nb3 +(x4-x1).*Nb4;
  ## Volume of the reference tetrahedron
  Kkvolume = 1/6;
  
  b(:,:) = Kkvolume * weight * detJ;
  
endfunction

function [b] = shg(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
  ## COMPUTE GRADIENT OF SHAPE FUNCTIONS
  

  Nb2 = y1.*(z3-z4) + y3.*(z4-z1) + y4.*(z1-z3);
  Nb3 = y1.*(z4-z2) + y2.*(z1-z4) + y4.*(z2-z1);
  Nb4 = y1.*(z2-z3) + y2.*(z3-z1) + y3.*(z1-z2);
  
  ## Determinant of the Jacobian of the 
  ## transformation from the base tetrahedron
  ## to the tetrahedron K  
  detJ = (x2-x1).*Nb2 +(x3-x1).*Nb3 +(x4-x1).*Nb4;
 
  ## Shape function gradients follow
  ## First index represents space direction
  ## Second index represents the shape function
  ## Third index represents the tetrahedron number
  b(1,1,:) = (y2.*(z4-z3) + y3.*(z2-z4) + y4.*(z3-z2))./ detJ;
  b(2,1,:) = (x2.*(z3-z4) + x3.*(z4-z2) + x4.*(z2-z3))./ detJ;
  b(3,1,:) = (x2.*(y4-y3) + x3.*(y2-y4) + x4.*(y3-y2))./ detJ;
  
  b(1,2,:) = ( Nb2 ) ./ detJ;
  b(2,2,:) = (x1.*(z4-z3) + x3.*(z1-z4) + x4.*(z3-z1)) ./ detJ;
  b(3,2,:) = (x1.*(y3-y4) + x3.*(y4-y1) + x4.*(y1-y3)) ./ detJ;
  
  b(1,3,:) = ( Nb3 ) ./ detJ;
  b(2,3,:) = (x1.*(z2-z4) + x2.*(z4-z1) + x4.*(z1-z2)) ./ detJ;
  b(3,3,:) = (x1.*(y4-y2) + x2.*(y1-y4) + x4.*(y2-y1)) ./ detJ;
  
  b(1,4,:) = ( Nb4) ./ detJ;
  b(2,4,:) = (x1.*(z3-z2) + x2.*(z1-z3) + x3.*(z2-z1)) ./ detJ;
  b(3,4,:) = (x1.*(y2-y3) + x2.*(y3-y1) + x3.*(y1-y2)) ./ detJ;
endfunction

%!shared mesh,wjacdet,shg,shp
% x = y = z = linspace(0,1,2);
% [mesh] = MSH3Mstructmesh(x,y,z,1,1:6)
% [wjacdet] =  MSH3Mgeomprop(mesh,"wjacdet")
% [shg] =  MSH3Mgeomprop(mesh,"shg")
% [shp] =  MSH3Mgeomprop(mesh,"shp")
%!test
% assert(columns(mesh.t),columns(wjacdet))
%!test
% assert(size(shg),[3 4 6])
%!test
% assert(shp,eye(4))
%!test
% fail(MSH3Mgeomprop(mesh,"samanafattababbudoiu"),"warning","Unexpected value in passed string. Empty vector passed as output.")