################################################################################
# File         : numerical_phase_velocity.jl
# Author       : Sandeep Koranne (C) 2017, All rights reserved.
# Purpose      : HW Assignment code for calculating numerical phase velocity
#              : using Newton's method
################################################################################
# Newton method for root
# x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
# question: how to get f'(x_n), Soln: use Finite Difference Approximation
# f'(x) = [f(x+h)-f(x-h)]/2h centered difference
function finite_difference_derivative(f,x,h)
    return ( ( f(x+h) - f(x-h) )/ (2*h) );
end

# Helper functions
function A(h,theta)
    return (0.5 * h * cos( theta ) );
end
function B(h,theta)
    return (0.5 * h * sin( theta ) );
end
function C(nu,N_PPW)
    return ( (1/(nu^2)) * ( sin( (pi*nu)/(N_PPW*1.0) )^2 ) );
end

function F(x,h,theta,nu,N_PPW)
    return ( ( sin( A(h, theta) * x)^2 ) + ( sin( B(h,theta) * x )^2 ) - C( nu, N_PPW ) );
end

function FD(x,h,theta,nu,N_PPW)
    return ( A(h,theta)*sin(2*A(h,theta)*x) + B(h,theta)*sin(2*B(h,theta)*x) );
end

# Or the explicit function derivative may be given

function NewtonRaphsonSolve( x0,theta,N_PPW,tolerance )
   nu = 0.5; # fixed for this problem
   h  = 1.0/N_PPW;
   MAX_ITERATION_COUNT = 1000;
   for i in 1:MAX_ITERATION_COUNT
       xn = x0 - F(x0,h,theta,nu,N_PPW)/FD(x0,h,theta,nu,N_PPW);
       #print("xn = "); println(xn);
       #print("x0 = "); println(x0);
       #print("f(x0) = "); println(F(x0,h,theta,nu,N_PPW));
       if( abs( xn-x0 ) < tolerance ) 
           #print("Solution found in : ");
           #println(i);
           return x0;
       end
       x0 = xn;
   end
end

function CalculateForAllAngles(N,N_PPW)
    THETA_MATRIX=zeros(N+1,2)
    for i in 1:N+1
        theta = (i-1)*2*pi/(N);
        THETA_MATRIX[i,1] = theta*360/(2*pi);
        THETA_MATRIX[i,2] = (2*pi)/NewtonRaphsonSolve( 2*pi, theta, N_PPW, 1e-6);
    end
    writedlm("N.txt", THETA_MATRIX, " ");
end

#sx = NewtonRaphsonSolve( 2*pi,0,5,1e-6);
#println("Newton Raphson solve: ");
#println( sx );
CalculateForAllAngles( 60, 20 );

