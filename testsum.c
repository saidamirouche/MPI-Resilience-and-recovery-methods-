#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <mpi-ext.h>
#define  ARRAYSIZE	16000000
#define  MASTER		0

float  data[ARRAYSIZE];
int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm);
void mpi_error_handler(MPI_Comm *comm, int *error_code,...);


int main (int argc, char *argv[])
{
  int   numtasks, taskid, rc, dest, offset, i, j, tag1,
    tag2, source, chunksize; 
  float mysum, sum;
  float update(int myoffset, int chunk, int myid);
  MPI_Status status;
  MPI_Comm world; /* a world comm for the work, w/o the spares */
  
  
  MPI_Init(&argc, &argv);
  rc=1;
  MPI_Errhandler new_eh;
 

  MPI_Comm_dup( MPI_COMM_WORLD, &world );

  MPI_Comm_create_errhandler(mpi_error_handler, &new_eh);
  
  MPI_Comm_set_errhandler(world, new_eh);
  //MPI_Comm_set_errhandler(MPI_COMM_WORLD, new_eh);
  MPI_Comm_size(world, &numtasks);
  
  MPI_Comm_rank(world,&taskid);
  printf ("MPI task %d has started...\n", taskid);

  
  MPI_Comm_size( world, &numtasks );
  MPI_Comm_rank( world, &taskid);
  chunksize = (ARRAYSIZE / numtasks);
  tag2 = 1;
  tag1 = 2;

  
  int victim;
  srand (time(NULL));
  victim = (rand()%(numtasks-1))+ 1;
  //victim = rand()%(numtasks-1);
  if( taskid == victim ) {
    printf( "Rank %04d: committing suicide\n", taskid );
    raise( SIGKILL );
    while(1); /* wait for the signal */
  }
    
  if (taskid == MASTER){

    sum = 0;
    for(i=0; i<ARRAYSIZE; i++) {
      data[i] =  i * 1.0;
      sum = sum + data[i];
    }

    printf("\n\n");
    offset = chunksize;
    for (dest=1; dest<numtasks; dest++) {
       MPI_Send(&offset, 1, MPI_INT, dest, tag1,world);
      /*  if( MPI_SUCCESS != rc ) {
	MPI_Comm_call_errhandler(world,rc);
      
	} */
       rc =  MPI_Send(&data[offset], chunksize, MPI_FLOAT, dest, tag2, world);
        if( MPI_SUCCESS != rc ) {
	  //MPI_Comm_call_errhandler(world,rc);
	} 
      printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset);
      offset = offset + chunksize;
    }

    offset = 0;
    mysum = update(offset, chunksize, taskid);

    for (i=1; i<numtasks; i++) {
      source = i;
       MPI_Recv(&offset, 1, MPI_INT, source, tag1, world, &status);
       /*   if( MPI_SUCCESS != rc ) {
	MPI_Comm_call_errhandler(world,rc);
      
	} */ 
       MPI_Recv(&data[offset], chunksize, MPI_FLOAT, source, tag2,
		     world, &status);
       /* if( MPI_SUCCESS != rc ) {
	MPI_Comm_call_errhandler(world,rc);
      
	} */ 
    }

    MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, MASTER, world);
    printf("Sample results: \n");
    offset = 0;
    for (i=0; i<numtasks; i++) {
      for (j=0; j<5; j++) 
	printf("  %e",data[offset+j]);
      printf("\n");
      offset = offset + chunksize;
    }
    printf("*** Final sum= %e ***\n",sum);

  }  
 

  if (taskid > MASTER) {

   
    source = MASTER;
    MPI_Recv(&offset, 1, MPI_INT, source, tag1, world, &status);
    MPI_Recv(&data[offset], chunksize, MPI_FLOAT, source, tag2, 
	     world, &status);

    mysum = update(offset, chunksize, taskid);

    dest = MASTER;
    MPI_Send(&offset, 1, MPI_INT, dest, tag1, world);
    MPI_Send(&data[offset], chunksize, MPI_FLOAT, MASTER, tag2, world);

    MPI_Reduce(&mysum, &sum, 1, MPI_FLOAT, MPI_SUM, MASTER, world);

  } 

  MPI_Comm_free( &world );

  MPI_Finalize();
  return EXIT_SUCCESS;
}


float update(int myoffset, int chunk, int myid) {
  int i; 
  float mysum;
  mysum = 0;
  for(i=myoffset; i < myoffset + chunk; i++) {
    data[i] = data[i] + i * 1.0;
    mysum = mysum + data[i];
  }
  printf("Task %d mysum = %e\n",myid,mysum);
  return(mysum);
}

void mpi_error_handler(MPI_Comm *comm, int *error_code, ...)
{
    
  int num_failures,i,mpi_mcw_rank, mpi_mcw_size;
  int loc_size,size,rank,rank_key;
  MPI_Group f_group,group_orig;
  MPI_Comm new_comm,com_fault_tolerant,everyone;
  int* failed_ranks;
  int* all_ranks;
  int myrank;
  int err[3];
  MPI_Comm_size(*comm, &loc_size);
  MPI_Comm_rank(*comm, &myrank);
  int rc;
    
  OMPI_Comm_failure_ack(*comm);
  OMPI_Comm_failure_get_acked(*comm, &f_group);
  MPI_Comm_group(*comm,&group_orig);
  
  MPI_Group_size(f_group, &num_failures);

  printf("nombre des echecs est %d \n", num_failures);
  failed_ranks = (int *)malloc(sizeof(int)*num_failures);
  all_ranks    = (int *)malloc(sizeof(int)*num_failures);

  for( i = 0; i < num_failures; ++i )
    {
      failed_ranks[i] = i;
    }
  MPI_Group_translate_ranks(f_group, num_failures, failed_ranks,group_orig, all_ranks);
  

  //------------------------------------------ Communicator maintenance------------------------------------------
  if(*error_code == MPI_ERR_PROC_FAILED ||
     *error_code == MPI_ERR_REVOKED
     )
    {
      //------------------------------------------ revoke ------------------------------------------
      if(*error_code != MPI_ERR_REVOKED)
	{
	  OMPI_Comm_revoke(*comm); }
	
      MPI_Comm_size(*comm, &loc_size);
      //------------------------------------------ shrink ------------------------------------------
      OMPI_Comm_shrink(*comm, &new_comm);

      com_fault_tolerant = new_comm;
      new_comm = MPI_COMM_NULL;
      MPI_Comm_size(com_fault_tolerant, &mpi_mcw_size);
      MPI_Comm_rank(com_fault_tolerant, &mpi_mcw_rank);

      printf("Processus %d en Global et %d en Shrinked comm \n", myrank,mpi_mcw_rank);

      rc = MPI_Comm_spawn("res",MPI_ARGV_NULL, 1, MPI_INFO_NULL, 0, com_fault_tolerant, &everyone, err );

      if (MPI_SUCCESS == rc) {
	
	fprintf(stdout,"done spawn \n");
	fflush(stdout);
      } 
      
      rc = MPI_Intercomm_merge(everyone, 1, &com_fault_tolerant);
      
      if (MPI_SUCCESS == rc) {
	
	fprintf(stdout,"done merge \n");
	fflush(stdout);
      } else {printf("no merge \n");
      }
      
      MPI_Comm_size(com_fault_tolerant, &size);
      
      printf("le taille est de %d  \n ",size);
      
      MPI_Comm_rank(com_fault_tolerant, &rank);
      
      //------------------------------------------remet l'ordre dans le nouveau comm ------------------------------------------
      if(rank<all_ranks[0])
	rank_key=rank;
      else
	rank_key=rank+1;
      MPI_Comm_split(com_fault_tolerant,0,rank_key,&new_comm);
      com_fault_tolerant=new_comm; 
      
      printf("c bon ");
    }

  fflush(NULL);

  return;
    
}

