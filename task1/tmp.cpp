#include <mpi.h>
#include <iostream>
#include <ctime>
#include <cmath>
using namespace std;  

int bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == root) {
        for (int i = 0; i <size; i++) {
            if (i != rank) {
                MPI_Send(data, count, datatype, i, 0, communicator);
            }
        }
    } else {
        MPI_Recv(data, count, datatype, root, 0, communicator,
            MPI_STATUS_IGNORE);
    }
    return 0;
}

int main(int argc, char **argv)
{
    int n=10,m=9, print=1;
    int *arr;

    int *a = new int[4];
    
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank==1){
        for(int i=0;i<4;i++)
            a[i]=i;
    }

    bcast(a,4,MPI_INT,1,MPI_COMM_WORLD);
    a[2]+=rank-2;
    cout<<endl<<a[2];
    int numRowToSend = (n/size);
    int realNumRow;

    MPI_Finalize();

    return 0;
}