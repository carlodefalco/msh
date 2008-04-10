function [msh] = MSH2Mequalizemesh(msh)


  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{Ax}, @var{Ay}, @var{bx}, @
  ## @var{by}]} = MSH2Mequalizemesh(@var{msh})
  ##
  ## To equalize the size of  triangle edges, apply a baricentric@
  ## regularization, i.e. move each node to the @
  ## center of mass of the patch of triangles to which it belongs.
  ##
  ## May be useful when distorting a mesh, 
  ## type @code{ demo MSH2Mequalizemesh } to see some examples. 
  ##
  ## @seealso{MSH2Mdisplacementsmoothing}
  ##
  ## @end deftypefn
  
  ## This file is part of 
  ##
  ##                   MSH - Meshing Software Package for Octave
  ##      -------------------------------------------------------------------
  ##              Copyright (C) 2007  Carlo de Falco
  ##              Copyright (C) 2007  Culpo Massimiliano
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

  x  = msh.p(1,:)';
  y  = msh.p(2,:)';

  dnodes = unique(msh.e(1:2,:)(:));
  varnodes = setdiff([1:columns(msh.p)],dnodes);

  Ax = spalloc(length(x),length(x),1);
  Ay = spalloc(length(x),length(x),1);

  ax = zeros(3,3,nel);
  ay = zeros(3,3,nel);

  for inode=1:3
    giinode(inode,:)=msh.t(inode,:);
    for jnode=1:3
      ginode(inode,jnode,:)=msh.t(inode,:);
      gjnode(inode,jnode,:)=msh.t(jnode,:);
    end
  end

  for ii=1:3  
    
    for jj=ii+1:3
      
      ax(ii,jj,:) = ax(jj,ii,:) = -ones(1,1,nel);
      ay(ii,jj,:) = ay(jj,ii,:) = -ones(1,1,nel);
      
      ax(ii,ii,:) -= ax(ii,jj,:);
      ax(jj,jj,:) -= ax(ii,jj,:);
      ay(ii,ii,:) -= ay(ii,jj,:);
      ay(jj,jj,:) -= ay(ii,jj,:);
      
    endfor
  endfor

  Ax = sparse(ginode(:),gjnode(:),ax(:));
  Ay = sparse(ginode(:),gjnode(:),ay(:));

  x(varnodes) = Ax(varnodes,varnodes) \ (-Ax(varnodes,dnodes)*x(dnodes));
  y(varnodes) = Ay(varnodes,varnodes) \ (-Ay(varnodes,dnodes)*y(dnodes));
  msh.p(1,:) = x';
  msh.p(2,:) = y';

%!demo
%! ### equalize a structured mesh without moving boundary nodes
%! msh = MSH2Mstructmesh(linspace(0,1,10),linspace(0,1,10),1,1:4,"random");
%! dnodes = MSH2Mnodesonsides(msh,1:4);
%! varnodes = setdiff([1:columns(msh.p)],dnodes);
%! x = msh.p(1,:)';
%! y = msh.p(2,:)';
%! msh = MSH2Mequalizemesh(msh);
%! triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');
%! pause(.01)

%!demo
%! ### distort a mesh on a square equalizing at each step
%! msh = MSH2Mstructmesh(linspace(0,1,10),linspace(0,1,10),1,1:4,"random");
%! dnodes = MSH2Mnodesonsides(msh,1:4);
%! varnodes = setdiff([1:columns(msh.p)],dnodes);
%! x = msh.p(1,:)';
%! y = msh.p(2,:)';
%! dx = dy = zeros(columns(msh.p),1);
%! dytot = dxtot = -.7*sin(x(dnodes).*y(dnodes)*pi/2);
%! Nsteps = 30;
%! for ii=1:Nsteps
%!  dx(dnodes) = dxtot;
%!  dy(dnodes) = dytot;
%!  [Ax,Ay] = MSH2Mdisplacementsmoothing(msh,1);
%!  dx(varnodes) = Ax(varnodes,varnodes) \ ...
%!      (-Ax(varnodes,dnodes)*dx(dnodes));
%!  dy(varnodes) = Ay(varnodes,varnodes) \ ...
%!  (-Ay(varnodes,dnodes)*dy(dnodes));
%!  msh.p(1,:) += dx'/Nsteps;
%!  msh.p(2,:) += dy'/Nsteps;
%! triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)','r');
%!  pause(.5)
%! x = msh.p(1,:)';
%! y = msh.p(2,:)';
%! msh = MSH2Mequalizemesh(msh);
%! hold on;triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');hold off
%!  pause(.5)
%! endfor