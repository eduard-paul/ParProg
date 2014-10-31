//#include "mpi.h"
//#include <iostream>
//#include <ctime>
//#include <iomanip>
//#include <stdio.h>
//using namespace std;
//
//int main(int argc, char **argv)
//{
//    MPI_Init(&argc, &argv);
//
//    /* int n = atoi(argv [1]);    
//    int m = atoi(argv [2]);*/
//    int n=10, m=5;
//    // код с использованием MPI
//    int *Mat;
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (rank==0){
//        srand(time(0));
//        //создание матрицы
//        Mat = new int [n*m]; 
//        //заполнение
//        for (int i = 0; i < n; i++)
//            for (int j = 0; j < m; j++)
//                Mat[i*m+j] = (rand() % 10 + 1);
//        // вывод матрицы
//        for (int i = 0; i < (n); i++)
//        {
//            for (int j = 0; j < m; j++)
//                //cout << Mat[i*m+j] << "   ";
//                    printf("%d   ", Mat[i*m+j]);
//            printf("\n");
//        }    
//    }
//
//    int *buf;
//    int q = n/size;
//    int count;
//
//    if (rank==0){
//        for (int i=1; i<size-1; i++){
//            MPI_Send( &Mat[i*q*m], q*m, MPI_INT, i,  MPI_ANY_TAG, MPI_COMM_WORLD);
//        }
//        if (size!=1)MPI_Send( &Mat[(size-1)*q*m], n*m-(size-1)*q*m, MPI_INT, size-1, size-1, MPI_COMM_WORLD);
//        buf=Mat;
//        count=q*m; 
//    }
//    else{
//        MPI_Status s;
//        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &s); //параметры сообщения без чтения
//        MPI_Get_count(&s, MPI_CHAR, &count);
//
//        buf= new int [count];
//        MPI_Recv(buf, count, MPI_INT,
//            MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &s);
//    }
//    int Max=buf[0];
//    printf("Max = %d   ",Max);
//    for (int i=1; i<count; i++){
//        if (buf[i]>Max) 
//            Max=buf[i];
//    }
//    
//    int MaxMax;
//    MPI_Reduce(&Max, &MaxMax, 1, MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
//    if (rank==0){
//        printf("Max = %d   ",MaxMax);
//    }
//
//    delete []Mat;
//
//    MPI_Finalize();
//    getchar();
//    return 0;
//}