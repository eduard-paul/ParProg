//// Сумма значений по столбцам матрицы.
//
//#include <mpi.h>
//#include <iostream>
//#include <ctime>
//#include <cmath>
//using namespace std;  
//
//int main(int argc, char **argv)
//{
//    int n=10,m=9, print=1;
//    int *arr;
//
//    if(argc==3){
//        n=atoi(argv[1]);
//        m=atoi(argv[2]);
//    }
//
//    if(argc==4){
//        n=atoi(argv[1]);
//        m=atoi(argv[2]);
//        print=atoi(argv[3]);
//    }
//
//    MPI_Init(&argc, &argv);
//
//    int rank, size;
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//
//    int numRowToSend = (n/size);
//    int realNumRow;
//    if( rank==0 )
//    {
//        // Формирование матрицы
//        srand(time(0));
//        arr = new int[n*m];
//        for (int i=0;i<n*m;i++)
//            arr[i]=rand()%(n*m);
//
//        // Вывод матрицы
//        if(print==1){
//            for(int i=0;i<n;i++){
//                for(int j=0;j<m;j++)
//                    printf("%d\t",arr[i*m+j]);
//                printf("\n");
//            }
//        }
//
//        // Рассылка
//        for (int i=1; i<size-1;i++)
//            MPI_Send(&arr[i*numRowToSend*m], numRowToSend*m, MPI_INT, i, i, MPI_COMM_WORLD);
//        // Послыка хвоста
//        if (size!=1)
//            MPI_Send(&arr[(size-1)*numRowToSend*m], n*m-(size-1)*numRowToSend*m, MPI_INT, size-1, size-1, MPI_COMM_WORLD);
//        realNumRow = numRowToSend;
//    } 
//    else {
//        // Прием 
//        MPI_Status s;        
//        if (rank!=size-1){
//            arr = new int[numRowToSend*m];
//            MPI_Recv(arr, numRowToSend*m, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &s);
//            realNumRow = numRowToSend;
//        }
//        // Прием хвоста
//        else{
//            realNumRow = n-(size-1)*numRowToSend;
//            arr = new int[realNumRow*m];
//            MPI_Recv(arr, realNumRow*m, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &s);            
//        }
//    }
//
//    // Локальные подсчеты
//    int* result = new int [m];
//    for(int i=0;i<m;i++) result[i]=0;
//    for(int i=0;i<realNumRow;i++)
//        for (int j=0;j<m;j++)
//            result[j]+=arr[i*m+j];
//
//    // Глобальный подсчет
//    int *finalResult; 
//    if(rank==0){
//        finalResult = new int[m];
//        for(int i=0;i<m;i++) finalResult[i]=0;
//    }
//
//    MPI_Reduce(result, finalResult, m, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
//
//    MPI_Finalize();
//
//    printf("\n");
//    if(rank==0)
//        for(int j=0;j<m;j++)
//            printf("%d\t",finalResult[j]);
//
//    if(rank==0){
//        delete []finalResult;}
//    delete []arr;
//    delete []result;
//    return 0;
//}