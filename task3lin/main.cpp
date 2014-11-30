#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <ctime>

unsigned int uintMask(int i) {
    unsigned int m = 1<<i;
    return m;
}

__int64 doubleMask(int i) {
    __int64 m = (__int64)1<<i;
    return m;
}

void UINT_MSD_Radix_Sort(unsigned int *a, int left, int right, int numBit = 31){
    if(left>=right) return;
    int l=left, r = right;
    while(left<right) {
        while(0==(uintMask(numBit)&(a[left]))) 
            left++;
        while(0!=((a[right])&(uintMask(numBit)))) 
            right--;
        if(left<right) {
            unsigned int tmp = a[left];
            a[left] = a[right];
            a[right] = tmp;
        }
    }
    if(numBit==0) return;
    UINT_MSD_Radix_Sort(a,l,right,numBit-1);
    UINT_MSD_Radix_Sort(a,left,r,numBit-1);
}

void _double_MSD_Radix_Sort(__int64 *a, int left, int right, int q, int numBit = 62){
    if(left>=right) return;
    int l=left, r = right;
    while(left<right) {
        if(q == 1){
            while((doubleMask(numBit)!=(doubleMask(numBit)&(a[left])))&&(left<r)) 
                left++;
            while((doubleMask(numBit)==((a[right])&(doubleMask(numBit))))&&(l<right)) 
                right--;
            if(left<right) {
                __int64 tmp = a[left];
                a[left] = a[right];
                a[right] = tmp;
            }
        } else {
            while((doubleMask(numBit)==(doubleMask(numBit)&(a[left])))&&(left<r)) 
                left++;
            while((doubleMask(numBit)!=((a[right])&(doubleMask(numBit))))&&(l<right)) 
                right--;
            if(left<right) {
                __int64 tmp = a[left];
                a[left] = a[right];
                a[right] = tmp;
            }
        }
    }
    if(numBit==0) return;
    _double_MSD_Radix_Sort(a,l,right,q,numBit-1);
    _double_MSD_Radix_Sort(a,left,r,q,numBit-1);
}

void double_MSD_Radix_Sort(double *a_, int left, int right){
    __int64 *a = (__int64*)a_;
    int numBit = 63;
    if(left>=right) return;
    int l=left, r = right;
    while(left<right) {
        while((doubleMask(numBit)==((a[left])&doubleMask(numBit)))&&(left<r)) 
            left++;
        while((doubleMask(numBit)!=((a[right])&(doubleMask(numBit))))&&(l<right)) 
            right--;
        if(left<right) {
            __int64 tmp = a[left];
            a[left] = a[right];
            a[right] = tmp;
        }
    }
    _double_MSD_Radix_Sort(a,l,right,0);
    _double_MSD_Radix_Sort(a,left,r,1);
}

void merge(double *a1, double *a2, double *b, int n1, int n2) {
    int k1=0, k2=0, k=0;

    while((k1<n1)&&(k2<n2)) {
        if(a1[k1]<=a2[k2]){
            b[k++]=a1[k1++];
        } else {
            b[k++]=a2[k2++];
        }
    }
    while (k1<n1) b[k++]=a1[k1++];
    while (k2<n2) b[k++]=a2[k2++];
}

void merge(double *a, int _k1, int _k2, int _k3) {
    int k1=_k1, k2=_k2, k=0, kk=k2;
    double *b1 = new double[_k3-_k1+1];
    while((k1<kk)&&(k2<_k3+1)) {
        if(a[k1]<=a[k2]){
            b1[k++]=a[k1++];
        } else {
            b1[k++]=a[k2++];
        }
    }
    while (k1<kk) b1[k++]=a[k1++];
    while (k2<_k3+1) b1[k++]=a[k2++];
    for(int i = _k1;i<_k3+1;i++) a[i]=b1[i-_k1];
    delete[] b1;
}

void localParSort(double* a, int start, int end, int threadNum=1){
    omp_set_dynamic(0); 
    omp_set_num_threads(threadNum);
    int N = end-start+1;
#pragma omp parallel
    {
        int q = omp_get_thread_num();
        int numCount = N/omp_get_num_threads() + 1;
        double_MSD_Radix_Sort(a,__min(q*numCount+start,end),__min((1+q)*numCount-1+start,end));
#pragma omp barrier 
#pragma omp master 
        if (omp_get_num_threads()==2) {
            int k = N/omp_get_num_threads() + 1;
            merge(a,start,start+k,end);
        }
        if (omp_get_num_threads()==4) {
            int k2 = N/omp_get_num_threads() + 1, k3 = 2*(N/omp_get_num_threads()+1),
                k4 = 3*(N/omp_get_num_threads() + 1);
#pragma omp single 
            merge(a,start+0,start+k2,start+k3-1);
#pragma omp single 
            merge(a,start+k3,start+k4,end);
#pragma omp barrier 
#pragma omp master 
            merge(a,start+0,start+k3,end);
        }
    } 
}

double log2( double n )  
{  
    return log( n ) / log( 2 );  
}

double rand(double LO, double HI){
    return (LO + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI-LO))));
}

using namespace std;

int main(int argc, char **argv)
{
    srand(time(0));
    int rank, size, threadNum;
    double time1, minV=-1000, maxV=1000;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double *a=NULL, *b;
    int N = atoi(argv[1]);
    threadNum = atoi(argv[2]);
    int mode = atoi(argv[3]);
    //int N = 18;
    if(rank==0){
        a = new double[N];
        b = new double[N];
        for(int i=0;i<N;i++){
            a[i]=rand(minV,maxV);
            b[i]=a[i];
        }        
        if(mode==1){
            double time2=-omp_get_wtime();
            double_MSD_Radix_Sort(b,0,N-1);
            time2+=omp_get_wtime();
            std::cout<<"Lenear version: "<<time2<<std::endl;
            MPI_Finalize();
            return 0;
        }
    }
    int *sendcounts = new int[size],*displs = new int[size];
    for(int i=0; i< size-1;i++){
        sendcounts[i]=N/size+1;
        displs[i] = i*(N/size+1);
    }    
    sendcounts[size-1]=N-((size-1)*(int)(N/size+1));
    displs[size-1]=((size-1)*(int)(N/size+1));
    double *recvbuf = new double[sendcounts[rank]];

    if(rank==0) time1=-omp_get_wtime();

    MPI_Scatterv(a,sendcounts,displs,MPI_DOUBLE,recvbuf,sendcounts[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    int currCount = sendcounts[rank];
    localParSort(recvbuf,0,currCount-1,threadNum);
    //Send&Merge
    bool empty = false;
    unsigned int urank = rank;
    for(int i=0;i<(int)(log2(size)+0.001);i++){
        if(!empty){
            if(urank==(urank|uintMask(i))){                
                empty = true;

                MPI_Send(recvbuf,currCount,MPI_DOUBLE,(int)(urank^uintMask(i)),0,MPI_COMM_WORLD); 

                delete[] recvbuf;
            }
            else{
                MPI_Status status;
                int recvCount;
                double *tmpRecv = new double[currCount*2];

                MPI_Recv(tmpRecv,N,MPI_DOUBLE,(int)(urank^uintMask(i)),0,MPI_COMM_WORLD,&status);

                MPI_Get_count( &status, MPI_DOUBLE, &recvCount );
                double *res = new double[currCount+recvCount];

                merge(recvbuf,tmpRecv,res,currCount,recvCount);

                currCount += recvCount;
                delete[] recvbuf;
                delete[] tmpRecv;
                recvbuf = new double[currCount];
                for(int i=0;i<currCount;i++) recvbuf[i]=res[i];
                delete[] res;
            }
        }
    }
    if(rank==0){
        time1+=omp_get_wtime();
        std::cout<<"Parallel version: "<<time1<<std::endl;
        if(N<30){
            for(int i=0;i<N;i++){
                std::cout<<recvbuf[i]<<", ";
            }std::cout<<std::endl;
            for(int i=0;i<N;i++){
                std::cout<<b[i]<<", ";
            }
        }
    }
    MPI_Finalize();
    return 0;
}
