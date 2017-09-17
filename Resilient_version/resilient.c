#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <signal.h>
#include <time.h>
#include <mpi-ext.h>
#include <stdlib.h>


int rank=MPI_PROC_NULL, verbose=0; 
char** gargv;
char p = '1';

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

#define maxn 12
int skip_first_suicide = 0;

void static_data(int first, int last, int rank, int size,double xlocal[(12/4)+2][12]) {         // Recovery schemes for static Data

  int i, j;

  skip_first_suicide = 1;

  for (i=1; i<=maxn/size; i++) 
	for (j=0; j<maxn; j++) 
	    xlocal[i][j] = rank;
  if (rank == size-1)
    for (j=0; j<maxn; j++) {
	xlocal[first][j] = -1;
    }
  else if (rank == 0)
    for (j=0; j<maxn; j++) {
	xlocal[last][j] = -1;
    }
}

void save_dat(int count, double  local[(12/4)+2][12],int rank) {           //  Fonctions save_date et recovery_data pour faire le Checkpoint et le restart
    int i,j;                                                        
    FILE *file;
char s[10];
    sprintf(s,"%d" , rank); 
    file = fopen(s, "w");


    
    fprintf(file, "%d \n", count );
    for (i = 1; i <=(12/4); i++) {
        for (j = 0; j < 12; j++)
            {
                fprintf (file, "%f\n", local[i][j]);
            }
    }
    fclose(file);
  
}

int recovery_data(double local[(12/4)+2][12],int rank , int dead) {
  int i,j,recov_it;
  FILE *fp;

  if (rank==dead)
    skip_first_suicide = 1;

  char s[10];
  sprintf(s,"%d" , rank); 
  fp = fopen(s, "r");
    
  if( NULL == fp )
    return 1;

  fscanf(fp,"%d",&recov_it);
      
  for (i = 1; i <=(12/4); i++) {
    for (j = 0; j < 12; j++) {
      if (!fscanf(fp,"%lf",&local[i][j]))
	break;
	      
    }
  }

  
  
  fclose(fp);
  return recov_it;
}
 
void LI_data(int first, int last, double xlocal[(12/4)+2][12]) {
  int k=0,i,j;
  double diffnorm;
  double     xnew[(12/3)+2][12];

  skip_first_suicide = 1;
                                                              // Linear interpolation for the dead process 
  do {
    k = k+1;
    diffnorm = 0;
    for (i=first; i<=last; i++) 
      for (j=1; j<maxn-1; j++) {
	xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
		      xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
	diffnorm += (xnew[i][j] - xlocal[i][j]) * 
	  (xnew[i][j] - xlocal[i][j]);
      }

    for (i=first; i<=last; i++) 
      for (j=1; j<maxn-1; j++) 
	xlocal[i][j] = xnew[i][j];
    diffnorm = sqrt(diffnorm);
  } while (diffnorm > 1.0e-13);

}



void reset_strategie(int dead, int first, int last, int rank, int size,double xlocal[(12/4)+2][12]) {
  

    if (rank==dead)
      static_data(first, last,  rank, size,xlocal);    // Reset the static Data for the dead process
}


void ER_strategie(int dead, int rank, double xlocal[(12/4)+2][12],MPI_Comm comm) {
  
  int itcnt;

  //       skip_first_suicide = 1;
       
       itcnt = recovery_data(xlocal,rank,dead);          // Recovery data in ER strategie by loading data from file


}

void LI_strategie(int dead, int first, int last, int rank, int size,double xlocal[(12/4)+2][12],MPI_Comm world) {
  int j;
  MPI_Status status;

  if(dead ==0 ) {                      // for the case of the borders 
    
  
   

    if(rank ==dead+1)
      MPI_Send( xlocal[1], maxn, MPI_DOUBLE, dead, 1, 
		world );


    if (rank==dead) {

      static_data(first, last, rank, size, xlocal);               // Reset the static data for the dead process     

      MPI_Recv( xlocal[maxn/size+1], maxn, MPI_DOUBLE, rank + 1, 1, 
		world, &status );

      for (j=0; j<maxn; j++) 
	xlocal[0][j] = -1;
      
     
    

      LI_data(first,last, xlocal);
    }
  }
  else if(dead==size-1) {                    // for the case of the borders 

   

    if(rank ==dead-1)
 
MPI_Send( xlocal[maxn/size], maxn, MPI_DOUBLE, dead, 0, 
		world );

    if (rank==dead) {

      static_data(first, last, rank, size, xlocal);               // Reset the static data for the dead process     

       MPI_Recv( xlocal[0], maxn, MPI_DOUBLE, rank - 1, 0, 
		world, &status );
      
     

      for (j=0; j<maxn; j++) 
	xlocal[maxn/size+1][j] = -1;
      
     
    

      LI_data(first,last, xlocal);
    }
  }
  

    
  else {
    
    
    if(rank== dead-1)
      MPI_Send( xlocal[maxn/size], maxn, MPI_DOUBLE, dead, 0, 
		world );
    if(rank ==dead+1)
      MPI_Send( xlocal[1], maxn, MPI_DOUBLE, dead, 1, 
		world );
    
    if(rank==dead) {                                                        // Exchange the data with the neighboors  

      static_data(first, last, rank, size, xlocal);               // Reset the static data for the dead process     

      MPI_Recv( xlocal[maxn/size+1], maxn, MPI_DOUBLE, rank + 1, 1, 
		world, &status );
      MPI_Recv( xlocal[0], maxn, MPI_DOUBLE, rank - 1, 0, 
		world, &status );
      
      LI_data(first,last, xlocal);                                      // Linear interpolation
    
    }
    
  }
}


		  
/* This example handles a 12 x 12 mesh, on 4 processors only. */

int main( argc, argv )
int argc;
char **argv;
{
  int        rank, value, size, errcnt, toterr, i, j, itcnt,recov_it;
    int        i_first, i_last,rc;
    MPI_Status status;
    double     diffnorm, gdiffnorm,recov;
    double     xlocal[(12/4)+2][12];
    double     xrec[(12/4)+2][12];
    double     xnew[(12/3)+2][12];
    MPI_Comm world,rworld;
    FILE *fp;
    MPI_Request request;
    double t,ts;
    clock_t begin,end;
    double t1,t2;
    int flag;
    int c,jeton;
    int dead, Dead;
    int dead_loc=0;
    char s[10];


   
    jeton = 1;
    gargv = argv;
    
    MPI_Init( &argc, &argv );


    c = atoi(argv[1]);
    int x[size];
    int y[size];

     MPI_Comm_get_parent( &world );
     
    if( MPI_COMM_NULL == world  ) {
      // on initialise notre communicateur 
      MPI_Comm_dup( MPI_COMM_WORLD, &world );
      MPI_Comm_size( world, &size );
      MPI_Comm_rank( world, &rank );
      MPI_Comm_set_errhandler( world, MPI_ERRORS_RETURN );
      //printf("im the dead one %d \n",rank);

    } else {
      // si le communicateur existe déja et le processus est spawné , il reprend le travail
      MPIX_Comm_replace( MPI_COMM_NULL, &world );
      
      MPI_Comm_size( world, &size );
      MPI_Comm_rank( world, &rank );
      MPI_Comm_set_errhandler( world, MPI_ERRORS_RETURN );
      dead_loc = rank;

      fp = fopen("info.txt","r");
      fscanf(fp,"%d %lf %lf \n", &itcnt,&t,&gdiffnorm);
      fclose(fp);


      fp = fopen("dead.txt", "a");
      fprintf(fp, "%e %d %lf %d\n",gdiffnorm, itcnt-1,t,rank); ;
      fclose(fp);

     
      if (c == 1)
	reset_strategie(rank, i_first, i_last, rank, size,xlocal);                    // Choosing the recovery strategie as a input parametre 
      else if (c==2) 
       
	ER_strategie(rank, rank,  xlocal,world);
       
   
      goto reprise;
    }

    i_first = 1;
    i_last  = maxn/size;
    if (rank == 0)        i_first++;
    if (rank == size - 1) i_last--;

    if (rank==0){
      sprintf(s,"dead.txt"); 
      fp = fopen(s, "w");
      fclose(fp);
    }
    
    int victim;                                                                     //Un processus qui se suicide
    srand (time(NULL));
    victim = (rand()%(size-1))+ 1;
    victim =1;
    t = 0;
    ts = 0;
    
    t1 = MPI_Wtime();

    for (i=1; i<=maxn/size; i++) 
      for (j=0; j<maxn; j++) 
	xlocal[i][j] = rank;
    for (j=0; j<maxn; j++) {
      xlocal[i_first-1][j] = -1;
	xlocal[i_last+1][j] = -1;
    }

    int informer;
    int me;
    itcnt = 0;
    
    do {

      
      if (c == 2)
	save_dat(itcnt,xlocal,rank);

      if((rank == 1)  && (itcnt == 100)) {

	if (c!=0) {
	  	if( 0 == skip_first_suicide ) {
	  exit(-1);
	}
	skip_first_suicide = 0;
	}

      }
      
      if (rank < size - 1) 
	rc = MPI_Send( xlocal[maxn/size], maxn, MPI_DOUBLE, rank + 1, 0, 
		       world );
     

    
      
      /* Send down unless I'm at the bottom */
      if (rank > 0) 
	rc = MPI_Send( xlocal[1], maxn, MPI_DOUBLE, rank - 1, 1, 
		       world );
  	
      if (rank < size - 1)   
	rc = MPI_Recv( xlocal[maxn/size+1], maxn, MPI_DOUBLE, rank + 1, 1, 
		       world, &status );

      flag = (rc == MPI_SUCCESS);                                      //*ULFM* Detects process failure and informe all the others about the error 
       MPIX_Comm_agree(world, &flag);
      

      if (rank > 0)
	rc =MPI_Recv( xlocal[0], maxn, MPI_DOUBLE, rank - 1, 0, 
		      world, &status );
     
	  
   
      flag = (rc == MPI_SUCCESS);                                      //*ULFM* Detects process failure and informe all the others about the error 
       MPIX_Comm_agree(world, &flag);
      
	/* Compute new values (but not on boundary) */
	itcnt ++;
	diffnorm = 0.0;
	for (i=i_first; i<=i_last; i++) 
	  for (j=1; j<maxn-1; j++) {
	    xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
			  xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
	    diffnorm += (xnew[i][j] - xlocal[i][j]) * 
	      (xnew[i][j] - xlocal[i][j]);
	  }
	/* Only transfer the interior points */
	for (i=i_first; i<=i_last; i++) 
	  for (j=1; j<maxn-1; j++) 
	    xlocal[i][j] = xnew[i][j];
       	
	
	rc = MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_SUM,
			    world );

	
	if( (MPI_ERR_PROC_FAILED == rc) || (MPI_ERR_PENDING == rc) ) {
	  dead_loc=0;

	  sprintf(s,"info.txt"); 
	  fp = fopen(s, "w");
	  fprintf(fp, "%d %lf %e\n", itcnt,t,gdiffnorm); ;
	  fclose(fp);
      

    	  flag = (rc == MPI_SUCCESS);
	  //	  MPIX_Comm_agree(world, &flag);
	  if( !flag ) {
	    MPIX_Comm_replace( world, &rworld );                      // Call the communicator fixing function ( shrink, spawn and reconstruct the COMM 
	    MPI_Comm_free( &world );                                  // and get back to recovery (reprise
	    world = rworld; 
	    goto reprise;
	  }
	} 
       	gdiffnorm = sqrt( gdiffnorm );

	t2 = MPI_Wtime();

	t =  (t2-t1);
	

	if (rank == 0) printf("%e %d %lf\n",gdiffnorm, itcnt,t);

    } while (gdiffnorm > 1.0e-13 && itcnt < 10000);

        MPI_Finalize( );
    return 0;
    
 reprise:
    

     MPI_Allreduce(&dead_loc, &dead, 1, MPI_INT, MPI_SUM,  world );

     // printf("%d :  iteration : %d  \n", rank,itcnt);

    
    i_first = 1;
    i_last  = maxn/size;
    if (rank == 0)        i_first++;
    if (rank == size - 1) i_last--;

 
    if(c==3)
	LI_strategie(dead, i_first, i_last, rank, size,xlocal,world);
      

   
    // itcnt = itcnt-2;                                 // redoing a previous iteration after recovery (reprise)
    
    do {
                                                      /* Send up unless I'm at the top, then receive from below */
                                                      /* Note the use of xlocal[i] for &xlocal[i][0] */
 
      if (c==2)
	save_dat(itcnt,xlocal,rank);                  // in case the recovery strategie chosen is ER , each process save its data in a file in each iteration

      
      if (rank < size - 1) 
	rc = MPI_Send( xlocal[maxn/size], maxn, MPI_DOUBLE, rank + 1, 0, 
		       world );
     
      
      /* flag = (rc == MPI_SUCCESS);
    
	 MPIX_Comm_agree(world, &flag);*/
  

      /* Send down unless I'm at the bottom */
      if (rank > 0) 
	rc = MPI_Send( xlocal[1], maxn, MPI_DOUBLE, rank - 1, 1, 
		       world );

     
      
      flag = (rc == MPI_SUCCESS);
      /*   if(!flag)
	   MPIX_Comm_agree(world, &flag);*/
      
     
      
      if (rank < size - 1) 
	rc = MPI_Recv( xlocal[maxn/size+1], maxn, MPI_DOUBLE, rank + 1, 1, 
		       world, &status );
    
	
      
      
      flag = (rc == MPI_SUCCESS);
      /*    if(!flag)
	    MPIX_Comm_agree(world, &flag);*/
     
	   
       
       if (rank > 0)
	rc = MPI_Recv( xlocal[0], maxn, MPI_DOUBLE, rank - 1, 0,    // Communications between neighboor processes for exchanging data
		       world, &status );

       

       
       flag = (rc == MPI_SUCCESS);
       //   if(!flag)
	 
	 //  MPIX_Comm_agree(world, &flag);                                 // *ULFM* Detects process failure and informe all the others about the error 

      
      
        if ((rank == 3) && (itcnt == 200) ) {                    
	                                                          
	if( 0 == skip_first_suicide ) {
	  exit(-1);
	}
	skip_first_suicide = 0;

	}

	if ((rank == 0) && (itcnt == 300) ) {                    
	                                                          
	if( 0 == skip_first_suicide ) {
	  exit(-1);
	}
	skip_first_suicide = 0;

	}
	
       
 
       /* -------------------------------1-----------------------------------
      
	  /* Compute new values (but not on boundary) */
       itcnt ++;
       diffnorm = 0;
       for (i=i_first; i<=i_last; i++) 
	 for (j=1; j<maxn-1; j++) {
	   xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
			 xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
	   diffnorm += (xnew[i][j] - xlocal[i][j]) * 
	     (xnew[i][j] - xlocal[i][j]);
	 }
       /* Only transfer the interior points */
	
       for (i=i_first; i<=i_last; i++) 
	 for (j=1; j<maxn-1; j++) 
	   xlocal[i][j] = xnew[i][j];


       


      
       rc = MPI_Allreduce( &diffnorm, &gdiffnorm, 1, MPI_DOUBLE, MPI_SUM,
			   world );


       if( (MPI_ERR_PROC_FAILED == rc) || (MPI_ERR_PENDING == rc) ) {

	 dead_loc=0;

	 sprintf(s,"info.txt"); 
	 fp = fopen(s, "w");
	 fprintf(fp, "%d %lf %e\n", itcnt,t,gdiffnorm); 

	 fclose(fp);
      

	 	 flag = (rc == MPI_SUCCESS);                                // Informe everybody about the fault and ask for recover steps
		 //	 MPIX_Comm_agree(world, &flag);
	 

	 if( !flag ) {
	   MPIX_Comm_replace( world, &rworld );                      // Call the communicator fixing function ( shrink, spawn and reconstruct the COMM 
	   MPI_Comm_free( &world );                                  // and get back to recovery (reprise)
	   world = rworld; 
	   goto reprise;
	 }
	  
       }

       
  if ((rank == 4) && (itcnt == 400) ) {                    
	                                                          
	if( 0 == skip_first_suicide ) {
	  exit(-1);
	}
	skip_first_suicide = 0;

      }
       
  
       
       gdiffnorm = sqrt( gdiffnorm );                           // Calculate the residus 

       t2 = MPI_Wtime();

       t = t2-t1;

       if (rank == 7) printf("%e %d %lf\n",gdiffnorm, itcnt,t);

       // MPI_Barrier(world);

       //  } while (itcnt < 100000);
         } while (gdiffnorm > 1.0e-13 && itcnt < 10000);


}
