// МЕТОДЫ ПЕРЕДАЧИ СООБЩЕНИЙ ([1] Лекция 5)
// Передача от одного всем (broadcast)

#include <mpi.h>
#include <iostream>
#include <ctime>
#include <cmath>
using namespace std;  

int tree_bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){

    int rank, size, retCode, *table;
    retCode=MPI_Comm_rank(MPI_COMM_WORLD, &rank); if(retCode!=0) return retCode;
    retCode=MPI_Comm_size(MPI_COMM_WORLD, &size);if(retCode!=0) return retCode;

    table = new int[size]; table[0]=-1;
    int k=1,q=0;
    for(int i=1;i<size;i++){
        if (q==k){
            q=0; k*=2;
        }
        table[i]=q++;
    }

    if (rank == root) {
        for(int i=0;i<size;i++){
            if(table[i]==rank){
                retCode=MPI_Send(data, count, datatype, i, 0, communicator);
                if(retCode!=0) return retCode;
            }
        } 
    } else {
        retCode=MPI_Recv(data, count, datatype, MPI_ANY_SOURCE, 0, communicator, MPI_STATUS_IGNORE);
        if(retCode!=0) return retCode;
        for(int i=0;i<size;i++){
            if(table[i]==rank){
                retCode=MPI_Send(data, count, datatype, i, 0, communicator);
                if(retCode!=0) return retCode;
            }
        } 
    }
    return retCode;
}

int bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){

    int rank, size, retCode;
    retCode=MPI_Comm_rank(MPI_COMM_WORLD, &rank); if(retCode!=0) return retCode;
    retCode=MPI_Comm_size(MPI_COMM_WORLD, &size);if(retCode!=0) return retCode;
    
    if (rank == root) {
        for (int i = 0; i <size; i++) {
            if (i != rank) {
                retCode=MPI_Send(data, count, datatype, i, 0, communicator);
                if(retCode!=0) return retCode;
            }
        }
    } else {
        retCode=MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
        if(retCode!=0) return retCode;
    }
    return retCode;
}

int main(int argc, char **argv)
{
   
    int N = atoi(argv[1]); 
    int m = atoi(argv[2]); 
    int root = atoi(argv[3]); 
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *a = new int[N];
    if(rank==root){
        for(int i=0;i<N;i++){
            a[i]=i;
        }
    }
  
    double bcast_time = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    bcast_time -= MPI_Wtime();

    for(int i=0;i<m;i++){
        MPI_Bcast(a,N,MPI_INT,root,MPI_COMM_WORLD);}

    MPI_Barrier(MPI_COMM_WORLD);
    bcast_time += MPI_Wtime();


    double my_bcast_time = 0;
    MPI_Barrier(MPI_COMM_WORLD); 
    my_bcast_time -= MPI_Wtime();

    for(int i=0;i<m;i++){
        bcast(a,N,MPI_INT,root,MPI_COMM_WORLD);}

    MPI_Barrier(MPI_COMM_WORLD);
    my_bcast_time += MPI_Wtime();

     double my_tree_bcast_time = 0;
    MPI_Barrier(MPI_COMM_WORLD); 
    my_tree_bcast_time -= MPI_Wtime();

    for(int i=0;i<m;i++){
        tree_bcast(a,N,MPI_INT,root,MPI_COMM_WORLD);}

    MPI_Barrier(MPI_COMM_WORLD);
    my_tree_bcast_time += MPI_Wtime();

    if(rank==0) cout<<"My broadcast time: "<<my_bcast_time<<endl<<"My tree_broadcast time: "<<my_tree_bcast_time<<endl<<"MPI broadcast time: "<<bcast_time<<endl;

    MPI_Finalize();

    return 0;
}