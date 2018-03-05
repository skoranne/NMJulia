################################################################################
# File   : heateqn.jl
# Author : Sandeep Koranne (C) 2018
# Purpose: Heterogenous heat equation for Homogenization experiment
#
################################################################################

using Gaston;
set(terminal="x11");
#using DiffEqOperators;

function TIME_VARYING_LEFT_BC(t)
    return 1+sin(t/10*pi);
end
function TIME_VARYING_RIGHT_BC(t)
    return 9+cos(t/10*pi);
end
function TIME_VARYING_SOURCE(t,x)
    return x*(1+sin(t/10*pi));
end

function CHECK_TIME_VARYING_FIELDS(L,N)
    U = zeros(N);
    for i in collect(1:N)
        U[i] = TIME_VARYING_SOURCE(i,0.5*L);
    end
    writedlm("CT.txt", U, " ");
end


function HETEROGENEOUS_HEAT_EQN()
    L  = 1.0;
    NUM_INTERVALS  = 10;
    Δx = L/NUM_INTERVALS;
    N  = NUM_INTERVALS-1; # number of internal nodes
    END_TIME = 100;
    TIME_ITERATIONS = 100;
    Δt = END_TIME/TIME_ITERATIONS;
    U = zeros(N+2,1); # temperature distribution
    UP= zeros(N+2,1); # temperature distribution of previous time step
    LEFT_BC  = 0.0;
    RIGHT_BC = 10.0;

    for i in collect(1:N+2)
        UP[i] = i*(RIGHT_BC - LEFT_BC)/(N+2);
    end
    writedlm("ui.txt", UP, " ");
    KAPPA = zeros(N,1);
    for i in collect(1:N)
        KAPPA[i] = 2+sin(2*pi*i*Δx/L); # since kappa(x) is positive
        println("KAPPA[", i, "] = ", KAPPA[i] );
    end
    RHS = zeros(N,1);
    for i in collect(1:N)
        if( i == 1 )
            RHS[i] = -LEFT_BC/(Δx)^2;
        end
        if( i == N )
            RHS[i] = -RIGHT_BC/(Δx)^2;
        end
    end
    #println("RHS = ", RHS);
    A = 1/(Δx)^2*(Tridiagonal( ones(N-1), -2*ones(N), ones(N-1)));
    #println("A=",A);
    x = A \ RHS ;
    #println("SOLN = ", A*x);
    writedlm("u.txt", x, " ");

    # for the varying coefficient case, we need the KAPPA values interpolated
    KAPPA_DL = zeros(N-1);
    KAPPA_D  = zeros(N);
    KAPPA_DU = zeros(N-1);

    for i in collect(2:N)
        KAPPA_DL[i-1] = 0.5*(KAPPA[i]+KAPPA[i-1]);
    end
    for i in collect(1:N-1)
        KAPPA_DU[i] = 0.5*(KAPPA[i]+KAPPA[i+1]);
    end
    for i in collect(1:N)
        if( i == 1 )
            KAPPA_D[i] = -1*( KAPPA_DL[1] + 0.5*(KAPPA[1]+KAPPA[2]));
        elseif( i == N )
            KAPPA_D[i] = -1*( KAPPA_DU[N-1] + 0.5*(KAPPA[N-1]+KAPPA[N]));
        else
            KAPPA_D[i] = -1*( KAPPA_DL[i-1] + KAPPA_DU[i-1] );
        end
    end
    writedlm("K.txt", KAPPA, " ");

    B = 1/(Δx)^2*(Tridiagonal( KAPPA_DL, KAPPA_D, KAPPA_DU ) );
    for t in collect(1:TIME_ITERATIONS)
        for i in collect(1:N)
            if( i == 1 )
                RHS[i] = -TIME_VARYING_LEFT_BC(t*Δt)/(Δx)^2;
            elseif( i == N )
                RHS[i] = -TIME_VARYING_RIGHT_BC(t*Δt)/(Δx)^2;
            else
                RHS[i] = TIME_VARYING_SOURCE(t*Δt,i*Δx);
            end
        end
        xk = B \ RHS;
        for i in collect(1:N)
            U[i+1] = UP[i+1]+ Δt*xk[i];
        end
        U[1]   = RHS[1]*(Δx)^2;
        U[N+2] = RHS[N]*(Δx)^2;
        for i in collect(1:N+2)
            UP[i] = U[i];
        end
    end
    writedlm("uk.txt", U, " ");    
end # end of function


#PLOT_HANDLE = figure();
# kplot = plot(KAPPA,title="KAPPA(x) coefficient",xlabel="x");
#uplot = plot(U,title="u(x) temperature distribution");
#printfigure( handle=PLOT_HANDLE,term="png",outputfile="x.png");

α  = 0.23;

#HETEROGENEOUS_HEAT_EQN();
CHECK_TIME_VARYING_FIELDS(1.0,100);

