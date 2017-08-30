############################################################
## File      : yee.jl
## Author    : Sandeep Koranne (C) All rights reserved. 2017
## Purpose   : Implementation of 1D Maxwell's equation with
##           : periodic boundary conditions.
## Method    : Yee's staggered grid.
############################################################
function yee(M,delta)
    #print("Running YEE FDTD 1D algorithm for M = ");
    #println(M);
    AX      = 0;
    BX      = 2*pi;
    AT      = 0;
    BT      = 20*pi;
    EZ_H    = M;
    HZ_H    = EZ_H-1;

    ERROR = 0;
    DX    = (BX-AX) / EZ_H;
    DT    = DX; # magic time steps
    TOTAL_TIME_STEPS = (BT-AT)/DT;
    EZ = zeros( EZ_H, 1 );
    HZ = zeros( HZ_H, 1 );
    EXACT_EZ = zeros( EZ_H, 1 );
    # Initial conditions
    for j in collect(1:EZ_H)
        EZ[j] = sin.( DX * (j-1) ) + DT*0.5 + delta * cos.((M/2)*DX*(j-1));
    end
    EZ[1] = EZ[EZ_H] = 0;
    for j in collect(1:HZ_H)
        HZ[j] = -sin.( DX * (j-0.5) + 0.5*DT );
    end
    for i=[1:EZ_H]
        EXACT_EZ[i] = sin.( DX*(i-1)  );
    end

    INITIAL_MATRIX = zeros( EZ_H-1, 4 );
    for j in collect(1:EZ_H-1)
        INITIAL_MATRIX[j,1] = DX * (j-1); # position vector
        INITIAL_MATRIX[j,2] = EXACT_EZ[j]; # exact field
        INITIAL_MATRIX[j,3] = EZ[j]; # numerical field
        INITIAL_MATRIX[j,4] = HZ[j];
    end
    writedlm("N.txt", INITIAL_MATRIX, " ");

    DT_DX_FACTOR = DT/DX;
    for n in collect(0:TOTAL_TIME_STEPS)
        if( false && mod.( n, 100 ) == 0 )
            print("Computed ");print(n);print(" of ");println(TOTAL_TIME_STEPS);
        end
        #EZ[1] = EZ[EZ_H] = sin.( DT*n +0.5*DT );
        for i=[2:EZ_H-1]
            EZ[i] = EZ[i] - DT_DX_FACTOR*(HZ[i]-HZ[i-1]);
        end
        for i=[1:HZ_H]
            HZ[i] = HZ[i] - DT_DX_FACTOR*(EZ[i+1]-EZ[i]);
        end
        EZ[1] = EZ[EZ_H] = EZ[EZ_H] - DT_DX_FACTOR*(HZ[1]-HZ[HZ_H]);
        EZ[1] = EZ[EZ_H] = sin.( DT*n + 0.5*DT);
        ##############################
        # Debug utility
        #
        if( n == -1 )
            for i=[1:EZ_H]
                EXACT_EZ[i] = sin.( DX*(i-1) + n*DT-0.5*DT  );
            end
            CHECK_MATRIX = zeros( EZ_H-1, 4 );
            for j in collect(1:EZ_H-1)
                CHECK_MATRIX[j,1] = DX * (j-1); # position vector
                CHECK_MATRIX[j,2] = EXACT_EZ[j]; # exact field
                CHECK_MATRIX[j,3] = EZ[j]; # numerical field
                CHECK_MATRIX[j,4] = HZ[j];
            end
            writedlm("N.txt", CHECK_MATRIX, " ");
        end
    end
    for i=[1:EZ_H]
        EXACT_EZ[i] = sin.( DX*(i-1) + 20*pi-0.5*DT );
    end
    ERROR = 0;
    for i in collect(1:EZ_H)
        value = ( EXACT_EZ[i] - EZ[i] );
        ERROR += value^2;
    end
    ERROR = sqrt.(DX*ERROR);
    #print(EZ);
    #print(EXACT_EZ);
    SUMMARY_MATRIX = zeros( EZ_H-1, 4 );
    for j in collect(1:EZ_H-1)
        SUMMARY_MATRIX[j,1] = DX * (j-1); # position vector
        SUMMARY_MATRIX[j,2] = EXACT_EZ[j]; # exact field
        SUMMARY_MATRIX[j,3] = EZ[j]; # numerical field
        SUMMARY_MATRIX[j,4] = HZ[j];
    end
    writedlm("M.txt", SUMMARY_MATRIX, " ");
    return ERROR;
end

for i in collect(3:9)
    M = 2^i;
    error = yee( M , 1e-3 );
    print("Error for M = ");
    print(M);
    print(" = " );
    println(error);
end
