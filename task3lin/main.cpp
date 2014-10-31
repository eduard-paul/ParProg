#include <iostream>
#include <omp.h>

unsigned int uintMask(int i) {
    unsigned int m = 1<<i;
    return m;
}

__int64 doubleMask(int i) {
    __int64 m = (__int64)1<<i;
    return m;
}

void UINT_LSD_Radix_Sort(unsigned int *a, int n){
    int numBit = 0;
    unsigned int *b1 = new unsigned int[n];
    unsigned int *b2 = new unsigned int[n];
    while(numBit < 32){
        int q1 = 0, q2 = 0;
        for(int i=0; i<n; i++)
            if(0==(uintMask(numBit)&(a[i]))) 
                b1[q1++]=a[i];
            else
                b2[q2++]=a[i];
        for(int i=0;i<q1;i++) a[i]=b1[i];
        for(int i=0;i<q2;i++) a[q1+i]=b2[i];
        numBit++;
    }
    delete []b1;
    delete []b2;
}

void double_LSD_Radix_Sort(double *a_, int left, int right){
    __int64 *a = (__int64*)a_;
    int numBit = 0;
    double *b1 = new double[right-left+1];
    double *b2 = new double[right-left+1];

    int q = 0, q2 = 0;
        for(int i=left; i<right+1; i++)
            if(doubleMask(63)==(doubleMask(63)&(a[i]))) 
                b1[q++]=a_[i];
            else
                b2[q2++]=a_[i];
        for(int i=0;i<q;i++) a_[i+left]=b1[i];
        for(int i=0;i<q2;i++) a_[q+i+left]=b2[i];

    while(numBit < 63){
        int q1 = 0, q2 = 0;
        for(int i=left; i<q+left; i++)
            if(doubleMask(numBit)==(doubleMask(numBit)&(a[i]))) 
                b1[q1++]=a_[i];
            else
                b2[q2++]=a_[i];

        for(int i=0;i<q1;i++) a_[i+left]=b1[i];
        for(int i=0;i<q2;i++) a_[q1+i+left]=b2[i];
        numBit++;
    }
    numBit=0;
    while(numBit < 63){
        int q1 = 0, q2 = 0;
        for(int i=q+left; i<right+1; i++)
            if(doubleMask(numBit)!=(doubleMask(numBit)&(a[i]))) 
                b1[q1++]=a_[i];
            else
                b2[q2++]=a_[i];

        for(int i=0;i<q1;i++) a_[q+i+left]=b1[i];
        for(int i=0;i<q2;i++) a_[q+q1+i+left]=b2[i];
        numBit++;
    }
    delete []b1;
    delete []b2;
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

void merge(double *a, double *b1, int k1, int k2, int k3) {
    int k=k1, kk=k2;
            while((k1<kk)&&(k2<k3+1)) {
                if(a[k1]<=a[k2]){
                    b1[k++]=a[k1++];
                } else {
                    b1[k++]=a[k2++];
                }
            }
            while (k1<kk) b1[k++]=a[k1++];
            while (k2<k3+1) b1[k++]=a[k2++];
}

int main(int argc, char **argv)
{
    int N = 1000000;
    double *a = new double[N];
    for(int i=0;i<N;i++){
        a[i]=N/2-i;
    }
    double *b = new double[N];
    for(int i=0;i<N;i++){
        b[i]=N/2-i;
    }

    double *b1 = new double[N];
    omp_set_dynamic(0); 
    omp_set_num_threads(2);
    double time1=-omp_get_wtime();
#pragma omp parallel
{
        int q = omp_get_thread_num();
        int numCount = N/omp_get_num_threads() + 1;
        double_MSD_Radix_Sort(a,__min(q*numCount,N-1),__min((1+q)*numCount-1,N-1));
#pragma omp barrier 
#pragma omp master 
        if (omp_get_num_threads()==2) {
            int k = N/omp_get_num_threads() + 1;
            merge(a,b1,0,k,N-1);
        }
        if (omp_get_num_threads()==4) {
            int k2 = N/omp_get_num_threads() + 1, k3 = 2*(N/omp_get_num_threads()+1),
                k4 = 3*(N/omp_get_num_threads() + 1);
#pragma omp single 
            merge(a,b1,0,k2,k3-1);
#pragma omp single 
            merge(a,b1,k3,k4,N-1);
#pragma omp barrier 
#pragma omp master 
            merge(b1,a,0,k3,N-1);
        }
}   
    time1+=omp_get_wtime();
    std::cout<<time1<<std::endl;

    double time2=-omp_get_wtime();
    double_MSD_Radix_Sort(b,0,N-1);
    time2+=omp_get_wtime();
    std::cout<<time2<<std::endl;

    //for(int i=0;i<N;i++){
    //    std::cout<<a[i]<<", ";
    //}std::cout<<std::endl;
    //for(int i=0;i<N;i++){
    //    std::cout<<b1[i]<<", ";
    //}
    system("pause");
    return 0;
}
