# This is a comment

# Euler integration for f'(x) = 3x^2, the solution should be x^3
# let f(0) = 0, f(1) = 1
# finite difference method is U_i - U_i-1 / h = f
# this gives U_i = f(x_i) * h + U_i-1 
# Or Implicit Backward Euler: U_{n+1} - U_{n} = f(x_n+1)*h
# this is in general non-linear which we have to solve for U_{n+1}
#f(x) = 3*x^2;
#TF(x) = x^3;
f(x) = -10*x;
TF(x) = exp(-10*x); 
N = 20;

# Newton method for root
# x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
# question: how to get f'(x_n), Soln: use Finite Difference Approximation
# f'(x) = [f(x+h)-f(x-h)]/2h centered difference
function finite_difference_derivative(f,x,h)
    return ( ( f(x+h) - f(x-h) )/ (2*h) );
end

function NewtonRaphsonSolve( f,x0,tolerance )
   MAX_ITERATION_COUNT = 1000;
   h = 0.001;
   for i in 1:MAX_ITERATION_COUNT
       xn = x0 - f(x0)/finite_difference_derivative( f,x0,h );
       print("xn = "); println(xn);
       print("x0 = "); println(x0);
       print("f(x0) = "); println(f(x0));
       if( abs( xn-x0 ) < tolerance ) 
           print("Solution found in : ");
           println(i);
           return x0;
       end
       x0 = xn;
   end
end

tbf(x) = 12-3*x^2;
sx = NewtonRaphsonSolve( tbf, 1.0, 0.01 );
println("Newton Raphson solve: ");
println( sx );


println("FD derivative of x^3 at x=2, h=0.01");
fdx = finite_difference_derivative( TF, 2, 0.01);
println(fdx);

function FE( A, B, DA, DB, N, f, TF )
    dh = 1/N;
    U = zeros( N, 1);
    U[1] = DA;
    U[N] = DB;
    for i in 2:N-1
        U[i] = f(U[i-1]) * dh + U[i-1];
    end
    TS = zeros( N,1);
    TS[1] = DA;
    TS[N] = DB;
    for i in 2:N-1
        TS[i] = TF( i*dh);
    end
    println("U");
    println(U);
    println("TS");
    println(TS);
    println("Norm");
    println(norm(U-TS));
end

function BE( A, B, DA, DB, N, f, TF )
    dh = 1/N;
    U = zeros( N, 1);
    U[1] = DA;
    U[N] = DB;
    # U_{n+1} - U_{n} = f(x_{n+1})*h
    for i in 2:N-1
        function fornewton(x) 
            return ( x-U[i-1] - dh*f(x) );
        end;
        U[i] = NewtonRaphsonSolve( fornewton, U[i-1], 1e-12 );
        #U[i] = f(i*dh) * dh + U[i-1];
    end
    TS = zeros( N,1);
    TS[1] = DA;
    TS[N] = DB;
    for i in 2:N-1
        TS[i] = TF( i*dh);
    end
    println("U");
    println(U);
    println("TS");
    println(TS);
    println("Norm");
    println(norm(U-TS));
end

FE( 0, 1, 1, exp(-10), 500, f, TF );
println("Now trying BE");
BE( 0, 1, 1, exp(-10), 500, f, TF );



