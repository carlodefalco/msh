function [Ax,Ay,Bx,By] = MSH2Mjigglemesh(msh, factor)


  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{Ax}, @var{Ay}, @var{bx}, @
  ## @var{by}]} @
  ## = MSH2Mjigglemesh(@var{msh}, @var{factor})
  ##
  ## To equalize the size of  triangle edges, set a spring of equilibrium
  ## length @var{factor}*@var{area} along each edge of the mesh and solve for
  ## equilibrium.
  ##
  ## This function builds matrices containing the resulting
  ## (linearized) equation for x and y coordinates of each mesh node.
  ## Boundary conditions enforcing the displacement (Dirichlet type
  ## problem) or the force (Neumann type) at the boundary must be added
  ## separately to make the system solvable.
  ##
  ## May be useful when distorting a mesh, e.g.:
  ##
  ## @example
  ##  
  ## msh = MSH2Mstructmesh(linspace(0,1,10),@
  ## linspace(0,1,10),@
  ## 1,1:4,"random");
  ## dnodes = MSH2Mnodesonsides(msh,1:4);
  ## varnodes = setdiff([1:columns(msh.p)],dnodes);
  ## xd = msh.p(1,dnodes)';
  ## yd = msh.p(2,dnodes)';
  ## dx = dy = zeros(columns(msh.p),1);
  ## dytot = dxtot = -.4*sin(xd.*yd*pi/2);
  ## Nsteps = 10;
  ## for ii=1:Nsteps
  ##  dx(dnodes) = dxtot;
  ##  dy(dnodes) = dytot;
  ##  [Ax,Ay] = MSH2Mdisplacementsmoothing(msh,1);
  ##  dx(varnodes) = Ax(varnodes,varnodes) \ ...
  ##      (-Ax(varnodes,dnodes)*dx(dnodes));
  ##  dy(varnodes) = Ay(varnodes,varnodes) \ ...
  ##  (-Ay(varnodes,dnodes)*dy(dnodes));
  ##  msh.p(1,:) += dx'/Nsteps;
  ##  msh.p(2,:) += dy'/Nsteps;
  ##  dx(dnodes) = 0;
  ##  dy(dnodes) = 0;
  ##  [AX, AY, BX,  BY] = MSH2Mjigglemesh(msh, .1);
  ##  dx(varnodes) = Ax(varnodes,varnodes) \ ...
  ##      (BX(varnodes)-Ax(varnodes,dnodes)*dx(dnodes));
  ##  dy(varnodes) = Ay(varnodes,varnodes) \ ...
  ##  (BY(varnodes)-Ay(varnodes,dnodes)*dy(dnodes));
  ##  msh.p(1,:) += dx';
  ##  msh.p(2,:) += dy';
  ##  triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');
  ##  pause(.01)
  ## endfor
  ## @end example
  ## @seealso{MSH2Mdisplacementsmoothing}
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
  ##   AUTHORS:
  ##   Carlo de Falco
  ##   Dublin City University
  ##   School of Mathemetical Sciences
  ##   Ireland
  ##
  ##   Culpo Massimiliano
  ##   Bergische Universitaet Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe fuer Angewandte MathematD-42119 Wuppertal  Gaussstr. 20 
  ##   D-42119 Wuppertal, Germany

  nel= columns(msh.t);

  x  = msh.p(1,:);
  y  = msh.p(2,:);

  dx  = (x(msh.t([1 2 3],:))-x(msh.t([2 3 1],:)));
  dy  = (y(msh.t([1 2 3],:))-y(msh.t([2 3 1],:)));
  dx2 = dx.^2;
  dy2 = dy.^2;

  l2  = dx2 + dy2;
  l   = sqrt(l2);
  area  = MSH2Mgeomprop(msh,'area');
  l0    = factor * sqrt(area) ([1 1 1],:);


  Ax = spalloc(length(x),length(x),1);
  Ay = spalloc(length(x),length(x),1);

  ax = zeros(3,3,nel);
  ay = zeros(3,3,nel);
  bx = zeros(3,nel);
  by = zeros(3,nel);

  for inode=1:3
    giinode(inode,:)=msh.t(inode,:);
    for jnode=1:3
      ginode(inode,jnode,:)=msh.t(inode,:);
      gjnode(inode,jnode,:)=msh.t(jnode,:);
    end
  end

  for ii=1:3  
    
    bx (ii,:) = ((l0(ii,:)-l(ii,:))./l(ii,:)) .* dx(ii,:);
    by (ii,:) = ((l0(ii,:)-l(ii,:))./l(ii,:)) .* dy(ii,:);

    for jj=ii+1:3
      
      ax(ii,jj,:) = ax(jj,ii,:) = reshape(1 +...
					  (l0(ii,:)./l(ii,:)).* ...
					  ((dx2(ii,:)./l2(ii,:)) -1 ) ,...
					  1,1,[]);
      
      ay(ii,jj,:) = ay(jj,ii,:) = reshape(1 +...
					  (l0(ii,:)./l(ii,:)).* ...
					  ((dy2(ii,:)./l2(ii,:)) -1 ) ,...
					  1,1,[]);
      
      ax(ii,ii,:) -= ax(ii,jj,:);
      ax(jj,jj,:) -= ax(ii,jj,:);
      ay(ii,ii,:) -= ay(ii,jj,:);
      ay(jj,jj,:) -= ay(ii,jj,:);
      
    endfor
  endfor

  Ax = sparse(ginode(:),gjnode(:),ax(:));
  Ay = sparse(ginode(:),gjnode(:),ay(:));
  Bx = sparse(giinode(:),1,bx(:));
  By = sparse(giinode(:),1,by(:));


endfunction

%!demo
