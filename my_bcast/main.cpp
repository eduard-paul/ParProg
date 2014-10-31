#include <mpi.h>
#include <stdio.h>
#include <ctime>
#include <iostream>
using namespace std;  

int brocast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator)
{
    int proc_rank, proc_kolvo, err;
    err=MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank); 
    if(err!=0) return err;
    err=MPI_Comm_size(MPI_COMM_WORLD, &proc_kolvo);
    if(err!=0) return err;

    if (proc_rank == root) 
    {
        for (int i = 0; i < proc_kolvo; i++) 
        {
            if (i != proc_rank) 
            {
                err=MPI_Send(data, count, datatype, i, 0, communicator);
                if(err!=0) return err;
            }
        }
    } 
    else 
    {
        err=MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
        if(err!=0) return err;
    }
    return err;
}

int main(int argc, char **argv)
{

    int n = atoi(argv[1]); 
    int cyclenum = atoi(argv[2]); 
    int root = atoi(argv[3]); 
    int proc_rank, proc_kolvo;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_kolvo);

    int *a = new int[n];
    if(proc_rank==root)
    {
        for(int i=0;i<n;i++)
        {    
            a[i]=i;
        }
    }

    double bcast_time;
    MPI_Barrier(MPI_COMM_WORLD);
    double bcast_time1 = MPI_Wtime();

    for(int i=0;i<cyclenum;i++){
        MPI_Bcast(a,n,MPI_INT,root,MPI_COMM_WORLD);}

    MPI_Barrier(MPI_COMM_WORLD);
    double bcast_time2 = MPI_Wtime();

    bcast_time=bcast_time2-bcast_time1;


    double brocast_time;
    MPI_Barrier(MPI_COMM_WORLD); 
    double brocast_time1 = MPI_Wtime();

    for(int i=0;i<cyclenum;i++) 
        brocast(a,n,MPI_INT,root,MPI_COMM_WORLD);
   
    MPI_Barrier(MPI_COMM_WORLD);
    double brocast_time2 = MPI_Wtime();

    brocast_time=brocast_time2-brocast_time1;

    if(proc_rank==0)
    {
        printf("\nmy broadcast=%f",brocast_time);
        printf("\nmpi broadcast=%f",bcast_time);
    }


    MPI_Finalize();

    return 0;
}