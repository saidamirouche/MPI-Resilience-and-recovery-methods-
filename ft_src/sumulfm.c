#include <mpi.h>
#include <mpi-ext.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>



int rank=MPI_PROC_NULL, verbose=0; 
char** gargv;

int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm) {
    MPI_Comm icomm,
             scomm, 
             mcomm; 
    MPI_Group cgrp, sgrp, dgrp;
    int rc, flag, rflag, i, nc, ns, nd, crank, srank, drank;

redo:
    if( comm == MPI_COMM_NULL ) {      //le nouveau processus attend son assignment de la part de rank 0
        MPI_Comm_get_parent(&icomm);
        scomm = MPI_COMM_WORLD;
        MPI_Recv(&crank, 1, MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
        if( verbose ) {
            MPI_Comm_rank(scomm, &srank);
            printf("Spawnee %d: crank=%d\n", srank, crank);
        }
    }
    else {  //la procedure de réparation pour les processus survivants 
      
      MPIX_Comm_shrink(comm, &scomm);    // shrinké le communicateur qui contient les echecs
      MPI_Comm_size(scomm, &ns);
      MPI_Comm_size(comm, &nc);
      nd = nc-ns; /* connaitre le nombre des processus tués  */
        if( 0 == nd ) {
            /* dans le cas ou y'a aucun processus tué , on fait rien */
            MPI_Comm_free(&scomm);
            *newcomm = comm;
            return MPI_SUCCESS;
        }
        
        MPI_Comm_set_errhandler( scomm, MPI_ERRORS_RETURN );

        rc = MPI_Comm_spawn(gargv[0], &gargv[1], nd, MPI_INFO_NULL,
                            0, scomm, &icomm, MPI_ERRCODES_IGNORE);
        flag = (MPI_SUCCESS == rc);
        MPIX_Comm_agree(scomm, &flag);
        if( !flag ) {
            if( MPI_SUCCESS == rc ) {
                MPIX_Comm_revoke(icomm);
                MPI_Comm_free(&icomm);
            }
            MPI_Comm_free(&scomm);
            if( verbose ) fprintf(stderr, "%04d: comm_spawn failed, redo\n", rank);
            goto redo;
        }

        /* se rappeller de l'ordre ancien pour le réassigner dans le nouveau communicator */
        MPI_Comm_rank(comm, &crank);
        MPI_Comm_rank(scomm, &srank);
        printf("mon rank dans shrink %d   ",srank);

	
        if(0 == srank) {
           
            MPI_Comm_group(comm, &cgrp);
            MPI_Comm_group(scomm, &sgrp);
            MPI_Group_difference(cgrp, sgrp, &dgrp);
            for(i=0; i<nd; i++) {
                MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
		// envoyer le nouveau assignement (rank) a le processus spawné 
                MPI_Send(&drank, 1, MPI_INT, i, 1, icomm);
            }
            MPI_Group_free(&cgrp); MPI_Group_free(&sgrp); MPI_Group_free(&dgrp);
        }
    }

    
    rc = MPI_Intercomm_merge(icomm, 1, &mcomm);
    rflag = flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(scomm, &flag);
    if( MPI_COMM_WORLD != scomm ) MPI_Comm_free(&scomm);
    MPIX_Comm_agree(icomm, &rflag);
    MPI_Comm_free(&icomm);
    if( !(flag && rflag) ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free(&mcomm);
        }
        if( verbose ) fprintf(stderr, "%04d: Intercomm_merge failed, redo\n", rank);
        goto redo;
    }

    rc = MPI_Comm_split(mcomm, 1, crank, newcomm);
    //réordre  le communicateur selon l'orignale dans World
    flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(mcomm, &flag);
    MPI_Comm_free(&mcomm);
    if( !flag ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free( newcomm );
        }
        if( verbose ) fprintf(stderr, "%04d: comm_split failed, redo\n", rank);
        goto redo;
    }

    /* restaurer l'erreur handler  */
    if( MPI_COMM_NULL != comm ) {
        MPI_Errhandler errh;
        MPI_Comm_get_errhandler( comm, &errh );
        MPI_Comm_set_errhandler( *newcomm, errh );
    }

    return MPI_SUCCESS;
}

#define  ARRAYSIZE	16000000
#define COUNT 1024
float  data[ARRAYSIZE];
float update(int myoffset, int chunk, int myid, float * tab) {
  int i; 
  float mysum;
  mysum = 0;
  for(i=myoffset; i < myoffset + chunk; i++) {
    mysum = mysum + tab[i];
  }
   printf("Task %d mysum  = %e   ",myid,mysum);
  return(mysum);
}

int main(int argc, char *argv[]) {
  int rank, np,rc;
    MPI_Comm world,rworld;
    double array[COUNT];
    //MPI_Init(NULL, NULL);
    gargv = argv;
    MPI_Init( &argc, &argv );
    char estr[MPI_MAX_ERROR_STRING]=""; int strl; /* error messages */
    double sum, mysum;
    int chunksize,i,offset,source,dest,tag1,tag2;
    MPI_Request send_request, recv_request;
    MPI_Status status;
    
    MPI_Comm_get_parent( &world );
    if( MPI_COMM_NULL == world ) {
      // on initialise notre communicateur 
      MPI_Comm_dup( MPI_COMM_WORLD, &world );
      MPI_Comm_size( world, &np );
      MPI_Comm_rank( world, &rank );
      MPI_Comm_set_errhandler( world, MPI_ERRORS_RETURN );
    } else {
      // si le communicateur existe déja et le processus est spawné , il reprend le travail
       MPIX_Comm_replace( MPI_COMM_NULL, &world );
      
      MPI_Comm_size( world, &np );
      MPI_Comm_rank( world, &rank );
      MPI_Comm_set_errhandler( world, MPI_ERRORS_RETURN );
     goto reprise;
    }
    

    int victim;                                                                     //Un processus qui se suicide
    srand (time(NULL));
    victim = (rand()%(np-1))+ 1;
    //victim = np-1;
    if((rank == victim) && (rank != 0)) {
      
      printf("\n Processus %d a commité un suicide  \n",rank);
      exit(0);
    }
    
    printf("mon ancien rank est %d   ", rank);

  chunksize = (ARRAYSIZE / np);


    

  if (rank == 0){
      
    sum = 0;
    for(i=0; i<ARRAYSIZE; i++) {
      data[i] =  i;
      sum = sum + data[i];
    }
    // printf("sum initial %e ",sum);
      
    offset = chunksize;
    for (dest=1; dest<np; dest++) {
      MPI_Send(&offset, 1, MPI_INT, dest, tag1,world);
       MPI_Send(&data[offset], chunksize, MPI_FLOAT, dest, tag2, world);
     
      offset = offset + chunksize;
    }

    
    
    offset = 0;
    mysum = update(offset, chunksize, rank,data);

   
      
    MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, world);
  
     
    printf("*** La somme finale = %e ***\n",sum);

  }  
   

 if (rank > 0) {
      
      source = 0;
      MPI_Recv(&offset, 1, MPI_INT, source, tag1, world,&status );
      MPI_Recv(&data[offset], chunksize, MPI_FLOAT, source, tag2 ,world,&status );

      

      mysum = update(offset, chunksize, rank,data);

       
      dest = 0;
     
       MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, world);

    } 

 
 rc = MPI_Barrier(world);
 if( (MPI_ERR_PROC_FAILED == rc) || (MPI_ERR_PENDING == rc) ) {
	MPIX_Comm_replace( world, &rworld );                      // appeller la fonction de réparation et remplacer le communicateur érroné avec le nouveau 
	MPI_Comm_free( &world );                                  // et commencer la partie de recovery 
 	world = rworld;
	goto reprise;
	}

 
     
 reprise:

 
    printf("mon nouveau rank est %d   \n\n", rank);
    
    chunksize = (ARRAYSIZE / np);


     if (rank == 0){
      	printf("\n broadcasting de nouveau ! \n");

        offset = chunksize;
      for (dest=1; dest<np; dest++) {
	MPI_Send(&offset, 1, MPI_INT, dest, tag1,world);
        rc = MPI_Send(&data[offset], chunksize, MPI_FLOAT, dest, tag2, world);
	printf("Envoie de %d elements a le rank %d \n",chunksize,dest);
	offset = offset + chunksize;
      }
      offset = 0;
      mysum = update(offset, chunksize, rank,data);


       MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, world);
       offset = 0;
       
       printf("*** La somme finale apres  shrink = %e ***\n",sum);

    }  
    if (rank > 0) {
      
      source = 0;
      MPI_Recv(&offset, 1, MPI_INT, source, tag1, world,&status );
      MPI_Recv(&data[offset], chunksize, MPI_FLOAT, source, tag2 ,world,&status );

      mysum = update(offset, chunksize, rank,data);
      dest = 0;
     
         MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, 0, world);

    }
    
      

    MPI_Comm_free( &world );
    MPI_Finalize();
    return EXIT_SUCCESS;
    
}

