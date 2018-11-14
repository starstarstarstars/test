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
//KM算法求二分图的最大完备匹配
//复杂度O(m*m*n) 

int n;  
int weight[MAX][MAX];           //权重,这是一个完备图，也就是二分图的左边节点集合X和右边节点集合Y之间任意两个点之间都有连通
//KM算法求二分图的最大权完美匹配，要求二分图本身存在完备匹配，如果没有是不能求得的
//但是KM算法求二分图最大权匹配，该图不一定完备，那么我们需要补充点，补充点之间的边权为0    
int lx[MAX],ly[MAX];                //定点标号  
bool sx[MAX],sy[MAX];          //记录寻找增广路时点集x，y里的点是否搜索过  
int match[MAX];                       //match[i]记录y[i]与x[match[i]]相对应  

bool search_path(int u) {          //给x[u]找匹配,这个过程和匈牙利匹配是一样的  
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
        //如果求最小匹配，则要将边权取反   
        if(!max_weight){ 
                for(int i=0;i<n;i++)  
                        for(int j=0;j<n;j++)  
                                weight[i][j]=-weight[i][j];  
        }  

        //初始化顶标，lx[i]设置为max(weight[i][j] | j=0,..,n-1 ), ly[i]=0;  
        for(int i=0;i<n;i++){  
                ly[i]=0;  
                lx[i]=-0x7fffffff;  
                for(int j=0;j<n;j++)  
                        if(lx[i]<weight[i][j])  
                                lx[i]=weight[i][j];  
        }  

        memset(match,-1,sizeof(match));  
        //不断修改顶标，直到找到完备匹配或完美匹配  
        for(int u=0;u<n;u++){   //为x里的每一个点找匹配  
                while(1){  
                        memset(sx,0,sizeof(sx));  
                        memset(sy,0,sizeof(sy));  
                        //x[u]在相等子图找到了匹配,继续为下一个点找匹配  
                        if(search_path(u))
                                break;  
                        //如果在相等子图里没有找到匹配，就修改顶标，直到找到匹配为止  
                        //首先找到修改顶标时的增量inc, min(lx[i] + ly [i] - weight[i][j],inc);,lx[i]为搜索过的点，ly[i]是未搜索过的点,因为现在是要给u找匹配，所以只需要修改找的过程中搜索过的点，增加有可能对u有帮助的边  
                        int inc=0x7fffffff;  
                        for(int i=0;i<n;i++)  
                                if(sx[i])  
                                        for(int j=0;j<n;j++)  
                                                if(!sy[j]&&((lx[i] + ly [j] - weight[i][j] )<inc))  
                                                        inc = lx[i] + ly [j] - weight[i][j] ;  
                         //找到增量后修改顶标，因为sx[i]与sy[j]都为真，则必然符合lx[i] + ly [j] =weight[i][j]，然后将lx[i]减inc，ly[j]加inc不会改变等式，但是原来lx[i] + ly [j] ！=weight[i][j]即sx[i]与sy[j]最多一个为真，lx[i] + ly [j] 就会发生改变，从而符合等式，边也就加入到相等子图中  
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