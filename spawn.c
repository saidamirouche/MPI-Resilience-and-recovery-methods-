#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define NUM_SPAWNS 1
char** gargv;

int main( int argc, char *argv[] )
{
    int np = NUM_SPAWNS;
    int errcodes[NUM_SPAWNS];
    MPI_Comm parentcomm, intercomm;
    int c,size,rank;
    
    gargv = argv;

    MPI_Init( &argc, &argv );

    gargv[1] = "1";

   
    printf("argument %s\n", argv[1]);
    
    MPI_Comm_get_parent( &parentcomm );

    // if(rank==0)
    //  if (parentcomm == MPI_COMM_NULL)
    
     if (argv[1] == "1")
      {      MPI_Comm_spawn( gargv[0], &gargv[1], np, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &intercomm, errcodes );
      
        printf("I'm the parent.\n");
    }
    else
    {
          printf("argument %s\n", argv[1]);

        printf("I'm the spawned.\n");
    }
    
    fflush(stdout);
    MPI_Finalize();
    return 0;
}
