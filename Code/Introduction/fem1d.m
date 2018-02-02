# File   : fem1d.m
# Author : Sandeep Koranne, for FEM class by Dr. Bokil
# 1D FEM code in Matlab
function y = f(x)
  y = pi^2*sin(pi*x);
endfunction

function y = exact(x)
  y = sin(pi*x);
endfunction

function y = hat1( x, xL, xR )
  y = (xR - x)/(xR-xL);
endfunction

function y = hat2( x, xL, xR )
  y = (x - xL)/(xR-xL);
endfunction

# use Simpson's rule to integrate
function y = int_hat1( xL, xR )
  mid_pt = (xL+xR)/2;
  y = (xR-xL)/6.0*(f(xL)*hat1(xL,xL,xR) + 4*f(mid_pt)*hat1(mid_pt,xL,xR) + f(xR)*hat1(xR,xL,xR));
endfunction

function y = int_hat2( xL, xR )
  mid_pt = (xL+xR)/2;
  y = (xR-xL)/6.0*(f(xL)*hat2(xL,xL,xR) + 4*f(mid_pt)*hat2(mid_pt,xL,xR) + f(xR)*hat2(xR,xL,xR));
endfunction

function [EXACT,U]=FEM(N, bc_left, bc_right)
h=1/N;
mesh=0:h:1;
ele=zeros(N,2);
for j = 1:N
  ele(j,1) = j;
  ele(j,2) = j+1;
  #printf("Fixing Element %d from (%f,%f)\n", j,mesh(ele(j,1)), mesh(ele(j,2)));
endfor
ele
mesh
K = 1/h*([1 -1;-1 1]);

LF = zeros(N,2);
for j = 1:N
  xL = mesh( ele(j, 1 ) );
  xR = mesh( ele(j, 2 ) );
  #printf("Load vector from %f to %f\n", xL, xR);
  LF(j,1) = int_hat1( xL, xR );
  LF(j,2) = int_hat2( xL, xR );
endfor


# Global assembly
F = zeros( N+1, 1);
A = zeros(N+1,N+1);

for j = 1:N
  nL = ele(j,1);
  nR = ele(j,2);
  #printf("Assembling for Element %d between node %d,%d\n", j,nL,nR);
  A(nL,nL) += K(1,1);
  A(nL,nR) += K(1,2);
  A(nR,nL) += K(2,1);
  A(nR,nR) += K(2,2);
  F(nL)    += LF(j,1);
  F(nR)    += LF(j,2);
endfor

# Impose Dirichlet boundary conditions
for j = 1:N+1
  A(1,j) = A(N+1,j) = 0;
endfor
A(1,1) = A(N+1,N+1) = 1.0;
F(1) = bc_left;
F(N+1) = bc_right;
U=A\F;
EXACT = exact(mesh);
endfunction

N=2;
[EXACT,U]=FEM(N,0,0);
EXACT
U
h=1/N;
x=0:h:1;
handle = figure(1);
plot(x,EXACT,'r',x,U,'b');
title 'EXACT and FEM solution for -u''=f, N=8';
xlabel 'x', ylabel 'u(x)';
print(handle,'-depsc2','plot2.eps');
quit;

