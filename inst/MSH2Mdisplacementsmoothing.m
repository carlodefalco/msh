function [Ax,Ay] = MSH2Mdisplacementsmoothing(msh, k)


  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{Ax},@var{Ay}]} = MSH2Mdisplacementsmoothing(@var{msh},@var{k})
  ##
  ## To displace the boundary of a 2D mesh, set a spring with
  ## force/length constant @var{k} along each edge and enforce
  ## equilibrium. This function builds matrices containing the resulting
  ## (linearized) equation for x and y coordinates of each mesh node.
  ## Boundary conditions enforcing the displacement (Dirichlet type
  ## problem) or the force (Neumann type) at the boundary must be added
  ## to make the system solvable, e.g.:
  ##
  ## @example
  ## msh = MSH2Mstructmesh(linspace(0,1,10),
  ##                      linspace(0,1,10),
  ##                      1,1:4,"left");
  ## dnodes = MSH2Mnodesonsides(msh,1:4);
  ## varnodes = setdiff([1:columns(msh.p)],dnodes);
  ## xd = msh.p(1,dnodes)';
  ## yd = msh.p(2,dnodes)';
  ## dx = dy = zeros(columns(msh.p),1);
  ## dxtot = dytot = -.5*sin(xd.*yd*pi/2);
  ## Nsteps = 10;
  ## for ii=1:Nsteps
  ##  dx(dnodes) = dxtot;
  ##  dy(dnodes) = dytot;
  ##  [Ax,Ay] = MSH2Mdisplacementsmoothing(msh,1);
  ##  dx(varnodes) = Ax(varnodes,varnodes) \ ...
  ##      (-Ax(varnodes,dnodes)*dx(dnodes));
  ##  dy(varnodes) = Ay(varnodes,varnodes) \ ...
  ##      (-Ay(varnodes,dnodes)*dy(dnodes));
  ##  msh.p += [ dx'/Nsteps; dy'/Nsteps ] ;
  ##  triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');
  ##  pause(.01)
  ## endfor
  ## @end example
  ##
  ## @end deftypefn
  
  ## This file is part of 
  ##
  ##                   MSH - Meshing Software Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco and Culpo Massimiliano
  ## 
  ##   MSH is free software; you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation; either version 2 of the License, or
  ##   (at your option) any later version.
  ## 
  ##   MSH is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ## 
  ##   You should have received a copy of the GNU General Public License
  ##   along with MSH; if not, write to the Free Software
  ##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  ##   USA
  ##
  ##
  ##   MAIN AUTHORS:
  ##   Carlo de Falco
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany

  x  = msh.p(1,:);
  y  = msh.p(2,:);

  dx2 = (x(msh.t([1 2 3],:))-x(msh.t([2 3 1],:))).^2;
  dy2 = (y(msh.t([1 2 3],:))-y(msh.t([2 3 1],:))).^2;

  l2  = dx2 + dy2;


  Ax = spalloc(length(x),length(x),1);
  Ay = spalloc(length(x),length(x),1);

  %for iel=1:columns(msh.t)

    ax = zeros(3,3,columns(msh.t));
    ay = zeros(3,3,columns(msh.t));
	
    for inode=1:3
      for jnode=1:3
	ginode(inode,jnode,:)=msh.t(inode,:);
	gjnode(inode,jnode,:)=msh.t(jnode,:);
      end
    end

    for ii=1:3  
      for jj=ii+1:3
	
	ax(ii,jj,:) = ax(jj,ii,:) = reshape(-k * dx2(ii,:)./l2(ii,:),1,1,[]);	
	ay(ii,jj,:) = ay(jj,ii,:) = reshape(-k * dy2(ii,:)./l2(ii,:),1,1,[]);

	ax(ii,ii,:) -= ax(ii,jj,:);
	ax(jj,jj,:) -= ax(ii,jj,:);
	ay(ii,ii,:) -= ay(ii,jj,:);
	ay(jj,jj,:) -= ay(ii,jj,:);

      endfor
    endfor

    Ax = sparse(ginode(:),gjnode(:),ax(:));
    Ay = sparse(ginode(:),gjnode(:),ay(:));
	 
 %endfor

endfunction

%!demo
%! msh = MSH2Mstructmesh(linspace(0,1,10),
%!                      linspace(0,1,10),
%!                      1,1:4,"left");
%! dnodes = MSH2Mnodesonsides(msh,1:4);
%! varnodes = setdiff([1:columns(msh.p)],dnodes);
%!
%! xd = msh.p(1,dnodes)';
%! yd = msh.p(2,dnodes)';
%!
%! dy = zeros(columns(msh.p),1);
%! dx = dy;
%!
%! dxtot = -.5*sin(xd.*yd*pi/2);
%! dytot = -.5*sin(xd.*yd*pi/2);
%!
%! Nsteps = 10;
%! for ii=1:Nsteps
%!  ii
%!  dx(dnodes) = dxtot;
%!  dy(dnodes) = dytot;
%!
%!  [Ax,Ay] = MSH2Mdisplacementsmoothing(msh,1);
%!  
%!  dx(varnodes) = Ax(varnodes,varnodes) \ ...
%!      (-Ax(varnodes,dnodes)*dx(dnodes));
%!  dy(varnodes) = Ay(varnodes,varnodes) \ ...
%!      (-Ay(varnodes,dnodes)*dy(dnodes));
%!
%!  msh.p(1,:) += dx'/Nsteps;
%!  msh.p(2,:) += dy'/Nsteps;
%!
%!    if mod(ii,2)==0
%!      triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');
%!      pause(.01)
%!    end
%! end
