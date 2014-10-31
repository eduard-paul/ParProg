//// lab1.cpp: определяет точку входа для консольного приложения.
////
//
//#include <iostream>
//#include <stdlib.h>
//#include "mpi.h"
//#include <time.h>
//using namespace std;
//
//double sum (double *ms, int rm)
//{
//    double res = 0;
//    for (int i = 0; i < rm; i++)
//    {
//        res += ms[i];
//    }
//    return res;
//}
//
//
//
//int main(int argc, char* argv[])
//{
//
//    int nTasks, rank;
//    MPI_Init(&argc,&argv);
//    MPI_Comm_size(MPI_COMM_WORLD,&nTasks);// получаем кол-во процессов
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);//получаем ранг текущего процесса
//    if (argc > 2)//инициализируем данные
//    {
//        int n = atoi(argv[1]);
//        int m = atoi(argv[2]);
//        int size = 0;
//        double **mt;
//        double *work;
//
//        double time1, time2,dt_psl,dt_prl;
//        double res_prl = 0,res_psl = 0;
//        if(((n*m)%nTasks) != 0)
//        {
//            size = (((n*m)/nTasks)+1)*nTasks;
//        }
//        else
//            size = n*m;
//
//        int klv = size/nTasks;
//        double *buf = new double[klv];//буффер для входящего сообщения
//        
//        double tmp2;//cумма элем в буффере
//        if (rank == 0)
//        {
//            srand((unsigned int)time(NULL));
//            mt = new double*[n];
//            for (int i = 0; i < n; i++)
//            {
//                mt[i] = new double[m];
//                for (int j = 0; j < m; j++)
//                    mt[i][j] = (rand()%100) + ((double)(rand()%100)/100);
//            }
//            //Копируем данные
//            work = new double[size];
//            for (int i = 0; i < n; i++)
//            {
//                for ( int j = 0;j < m; j++)
//                    work[(i*m)+j] = mt[i][j];
//            }
//
//            for (int i = (n * m); i < size; i++)	
//            {	
//                work[i] = work[(m * n) - 1];	
//            }
//
//            time1 = MPI_Wtime();
//            res_psl = sum(work,size);
//            time2 = MPI_Wtime();
//            dt_psl = time2 - time1;
//            
//
//        }
//        time1 = MPI_Wtime();
//        MPI_Scatter(work, klv, MPI_DOUBLE, buf, klv, MPI_DOUBLE, 0, MPI_COMM_WORLD);//рссылаем процессам блоки из work длиной klv
//        
//        tmp2 = sum(buf,klv);
//
//        MPI_Reduce(&tmp2, &res_prl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//coбираем информацию со всех процессов и суммируем
//
//        time2 = MPI_Wtime();
//        dt_prl = time2 - time1;
//
//        if (rank == 0)
//        {
//            cout << "Summa(posl) = "<<res_psl<<endl;
//            cout << "Time = "<<dt_psl<<endl;
//            cout << "Summa (paral) = "<<res_prl<<endl;
//            cout << "Time = "<<dt_prl<<endl;
//        }
//
//    }
//    else
//        if (rank == 0)
//            cout << "Error!"<<endl;
//
//    MPI_Finalize();
//    return 0;
//}