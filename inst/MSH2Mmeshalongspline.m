
function msh2 = MSH2Mmeshalongspline(xc,yc,Nnx,Nny,sigma)

  ## -*- texinfo -*-
  ## @deftypefn {Function File} {[@var{msh}]} =
  ## MSH2Mmeshalongspline(@var{xc}, @var{yc}, @var{Nnx}, @var{Nny},
  ## @var{sigma}) 
  ##
  ## Generates a structured mesh in a thin layer of size @var{sigma}
  ## around a natural cubic spline with control points @var{xc},
  ## @var{yc}. 
  ## The mesh will have @var{Nnx} nodes in the direction along the
  ## spline and @var{Nny} in the normal direction. 
  ## Be aware that if @var{sigma} is not much smaller than the curvature
  ## of the line the resulting mesh may be invalid. 
  ## 
  ## @seealso{MSH2Mstructmesh}
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
  ##   along with MSH; If not, see <http://www.gnu.org/licenses/>.
  ##
  ##
  ##   MAIN AUTHORS:
  ##   Culpo Massimiliano
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gaußstr. 20 
  ##   D-42119 Wuppertal, Germany
  ##
  ##   Carlo de Falco
  ##   Bergische Universität Wuppertal
  ##   Fachbereich C - Mathematik und Naturwissenschaften
  ##   Arbeitsgruppe für Angewandte MathematD-42119 Wuppertal  Gauupper. 20 
  ##   D-42119 Wuppertal, Germany

  s    =  [0; cumsum(sqrt(diff(xc).^2+diff(yc).^2))];
  xsPP =  csape ( s, xc , 'second', [0 0]);
  ysPP =  csape ( s, yc , 'second', [0 0]);
  
  ss   = linspace(0,s(end),Nnx);
  xs   = ppval(xsPP,ss);
  ys   = ppval(ysPP,ss);

  
  dxsPP = fnder(xsPP,1);
  dysPP = fnder(ysPP,1);


  nx = -ppval(dysPP,ss)';
  ny =  ppval(dxsPP,ss)';

  nx = nx ./ sqrt(nx.^2+ny.^2);
  ny = ny ./ sqrt(nx.^2+ny.^2);

  msh2 = MSH2Mstructmesh ([1:Nnx],linspace(-1,1,Nny),1,1:4);
  
  jj = (msh2.p(1,:));
  p(1,:) = xs(jj) + sigma*nx(jj)' .* msh2.p(2,:);
  p(2,:) = ys(jj) + sigma*ny(jj)' .* msh2.p(2,:);

  msh2.p = p;