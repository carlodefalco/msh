function [msh] = MSH2Mjigglemesh(msh, steps)


  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[newmsh]} = MSH2Mjigglemesh(@var{msh}, @var{steps})
  ##
  ## To equalize the size of  triangle edges, set a spring of rest
  ## length @var{factor}*@var{area} along each edge of the mesh and
  ## solve for static equilibrium.
  ##
  ## The non-linear eqautions of the system obtained are soved via a
  ## non-linear Gass-Seidel method. @var{step} is the number of steps of
  ## the method to be applied.
  ##
  ## May be useful when distorting a mesh, type @code{demo
  ## MSH2Mjigglemesh} to see some examples.
  ##
  ## @seealso{MSH2Mdisplacementsmoothing, MSH2Mequalizemesh}
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
  nnodes = columns(msh.p);

  x  = msh.p(1,:)';
  y  = msh.p(2,:)';

  dnodes = unique(msh.e(1:2,:)(:));
  vnodes = setdiff(1:nnodes,dnodes);

  %% find node neighbours XXX FIXME: should this go
  %% into MSH2Mtopprop ?
  sides = MSH2Mtopprop(msh,'sides');
  for inode = 1:nnodes
    neig{inode} = (sides(:, sides(1,:) == inode | sides(2,:) == inode))(:);
    neig{inode} (neig{inode} == inode) = [];
  endfor

  for istep = 1:steps
    for inode =vnodes
      xx = x(neig{inode}) * ones(size(neig{inode}))';
      lx = abs ( xx - xx' )(:);
      mx = ( xx + xx'  )(:)/2; 
      x(inode) = sum(mx.*lx)/sum(lx);

      yy = y(neig{inode}) * ones(size(neig{inode}))';
      ly = abs ( yy - yy' )(:);
      my = (yy + yy')(:)/2; 
      y(inode) = sum(my.*ly)/sum(ly);
    endfor
  endfor
  
  msh.p = [x';y'];
  
%!demo
%! ### distort a mesh on a square equalizing at each step
%! msh = MSH2Mstructmesh(linspace(0,1,10),linspace(0,1,10),1,1:4,"right");
%! dnodes = MSH2Mnodesonsides(msh,1:4);
%! varnodes = setdiff([1:columns(msh.p)],dnodes);
%! x = msh.p(1,:)';
%! y = msh.p(2,:)';
%! dx = dy = zeros(columns(msh.p),1);
%! dytot = dxtot = -.4*sin(x(dnodes).*y(dnodes)*pi/2);
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
%! msh = MSH2Mjigglemesh(msh,10);
%! hold on;triplot(msh.t(1:3,:)',msh.p(1,:)',msh.p(2,:)');hold off
%!  pause(.5)
%! endfor