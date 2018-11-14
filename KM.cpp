#include "mex.h"
#include <iostream>  
#include <cstdio>  
#include <memory.h>  
#include <algorithm>   

// using namespace std;  

#define MAX 1000  

#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]

#define	YP_OUT	plhs[0]
//KM�㷨�����ͼ������걸ƥ��
//���Ӷ�O(m*m*n) 

int n;  
int weight[MAX][MAX];           //Ȩ��,����һ���걸ͼ��Ҳ���Ƕ���ͼ����߽ڵ㼯��X���ұ߽ڵ㼯��Y֮������������֮�䶼����ͨ
//KM�㷨�����ͼ�����Ȩ����ƥ�䣬Ҫ�����ͼ��������걸ƥ�䣬���û���ǲ�����õ�
//����KM�㷨�����ͼ���Ȩƥ�䣬��ͼ��һ���걸����ô������Ҫ����㣬�����֮��ı�ȨΪ0    
int lx[MAX],ly[MAX];                //������  
bool sx[MAX],sy[MAX];          //��¼Ѱ������·ʱ�㼯x��y��ĵ��Ƿ�������  
int match[MAX];                       //match[i]��¼y[i]��x[match[i]]���Ӧ  

bool search_path(int u) {          //��x[u]��ƥ��,������̺�������ƥ����һ����  
        sx[u]=true;  
        for(int v=0; v<n; v++){  
                if(!sy[v] &&lx[u]+ly[v] == weight[u][v]){  
                        sy[v]=true;  
                        if(match[v]==-1 || search_path(match[v])){  
                                match[v]=u;  
                                return true;  
                        }  
                }  
        }  
        return false;  
}  

int Kuhn_Munkras(bool max_weight){
        //�������Сƥ�䣬��Ҫ����Ȩȡ��   
        if(!max_weight){ 
                for(int i=0;i<n;i++)  
                        for(int j=0;j<n;j++)  
                                weight[i][j]=-weight[i][j];  
        }  

        //��ʼ�����꣬lx[i]����Ϊmax(weight[i][j] | j=0,..,n-1 ), ly[i]=0;  
        for(int i=0;i<n;i++){  
                ly[i]=0;  
                lx[i]=-0x7fffffff;  
                for(int j=0;j<n;j++)  
                        if(lx[i]<weight[i][j])  
                                lx[i]=weight[i][j];  
        }  

        memset(match,-1,sizeof(match));  
        //�����޸Ķ��ֱ꣬���ҵ��걸ƥ�������ƥ��  
        for(int u=0;u<n;u++){   //Ϊx���ÿһ������ƥ��  
                while(1){  
                        memset(sx,0,sizeof(sx));  
                        memset(sy,0,sizeof(sy));  
                        //x[u]�������ͼ�ҵ���ƥ��,����Ϊ��һ������ƥ��  
                        if(search_path(u))
                                break;  
                        //����������ͼ��û���ҵ�ƥ�䣬���޸Ķ��ֱ꣬���ҵ�ƥ��Ϊֹ  
                        //�����ҵ��޸Ķ���ʱ������inc, min(lx[i] + ly [i] - weight[i][j],inc);,lx[i]Ϊ�������ĵ㣬ly[i]��δ�������ĵ�,��Ϊ������Ҫ��u��ƥ�䣬����ֻ��Ҫ�޸��ҵĹ������������ĵ㣬�����п��ܶ�u�а����ı�  
                        int inc=0x7fffffff;  
                        for(int i=0;i<n;i++)  
                                if(sx[i])  
                                        for(int j=0;j<n;j++)  
                                                if(!sy[j]&&((lx[i] + ly [j] - weight[i][j] )<inc))  
                                                        inc = lx[i] + ly [j] - weight[i][j] ;  
                         //�ҵ��������޸Ķ��꣬��Ϊsx[i]��sy[j]��Ϊ�棬���Ȼ����lx[i] + ly [j] =weight[i][j]��Ȼ��lx[i]��inc��ly[j]��inc����ı��ʽ������ԭ��lx[i] + ly [j] ��=weight[i][j]��sx[i]��sy[j]���һ��Ϊ�棬lx[i] + ly [j] �ͻᷢ���ı䣬�Ӷ����ϵ�ʽ����Ҳ�ͼ��뵽�����ͼ��  
                        for(int i=0;i<n;i++){  
                                if(sx[i])   //  
                                        lx[i]-=inc;  
                                if(sy[i])  
                                        ly[i]+=inc;  
                        }  
                }  

        }  
        int sum=0;  
        for(int i=0;i<n;i++)  
                if(match[i]>=0)  
                        sum+=weight[match[i]][i];  

        if(!max_weight)  
                sum=-sum;  
        return sum;  
}  
static void getMatch(
		   double	pout[],
		   double	t_weight[],
 		   double	sz[],
           bool     flag
		   )
{
    int i,j,sum,cnt;
    n = sz[0];
    cnt = sz[0];
    for(i=0;i<sz[0];i++)
        for(j=0;j<sz[0];j++)
        {
             weight[i][j] = t_weight[i*cnt+j];
//             pout[i*cnt+j] = t_weight[i*cnt+j];
        }           
   sum = Kuhn_Munkras(flag);
   for(i=0;i<MAX;i++)
       pout[i] = match[i];
   pout[i] = sum;
    
    return;
}
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )    
{ 
    double *yp; 
    double *t,*y;  
    bool  flag;
    YP_OUT = mxCreateDoubleMatrix( 1, MAX+1, mxREAL);  
    yp = mxGetPr(YP_OUT);	
    t = mxGetPr(prhs[0]); 
    y = mxGetPr(prhs[1]);
    flag = *(mxGetPr(prhs[2]));
    getMatch(yp,t,y,flag); 
    return;   
}