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
## @deftypefn {Function File} {[@var{msh}]} =
## MSH2Mmeshalongspline(@var{xc}, @var{yc}, @var{Nnx}, @var{Nny},
## @var{sigma}) 
##
## Generates a structured mesh in a thin layer of size @var{sigma}
## sitting on a natural Catmull-Rom type cubic spline with control
## points @var{xc}, @var{yc}. 
## If @var{Nnx} and @var{Nny} are scalars,
## the mesh will have @var{Nnx} nodes in the direction along the
## spline and @var{Nny} in the normal direction. 
## If @var{Nnx} and @var{Nny} are vectors they will indicate
## the curvilinear coordinates of the mesh nodes. 
## Be aware that if @var{sigma} is not much smaller than the curvature
## of the line the resulting mesh may be invalid. 
## 
## @seealso{MSH2Mstructmesh}
## @end deftypefn

 
function msh2 = MSH2Mmeshalongspline(xc,yc,Nnx,Nny,sigma)

  if (nargin != 5)
    print_usage();
  else
    s    =  [0:length(xc)-1];
    xsPP =  catmullrom ( s, xc );
    ysPP =  catmullrom ( s, yc );

    if (length(Nnx)>1)
      ss = Nnx(:).';
    else
      ss   = linspace(0,s(end),Nnx);
    endif
    xs   = ppval(xsPP,ss);
    ys   = ppval(ysPP,ss);
  
    dxsPP = fnder(xsPP,1);
    dysPP = fnder(ysPP,1);


    nx = -ppval(dysPP,ss)';
    ny =  ppval(dxsPP,ss)';
    
    nx = nx ./ sqrt(nx.^2+ny.^2);
    ny = ny ./ sqrt(nx.^2+ny.^2);

    if (length(Nny)>1)
      ssy = Nny(:).';
    else
      ssy = linspace(0,1,Nny);
    endif

    msh2 = MSH2Mstructmesh ([1:length(ss)], ssy, 1, 1:4);
    
    jj = (msh2.p(1,:));
    p(1,:) = xs(jj) + sigma*nx(jj)' .* msh2.p(2,:);
    p(2,:) = ys(jj) + sigma*ny(jj)' .* msh2.p(2,:);
    
    msh2.p = p;
  endif

endfunction