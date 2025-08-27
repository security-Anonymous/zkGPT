#include "hyrax.hpp"
#include "timer.hpp"
#include <cmath>
#include <thread>
#include <vector>
#include <mutex>
#include <iostream>
#include <chrono>
using namespace std;
using namespace mcl::bn;

const int MAX_MSM_LEN=1e4;
const int COMM_OPT_MAX=65536; //don't optimize if larger than this
const int logmax=16;  /// max number=2^18-1
const int block_num=5;  //2^80


G1 perdersen_commit(G1* g,ll* f,int n,G1* W)
{
    G1 ret;
    ret.clear();
    //timer t(true);
    //t.start();
    
    bool *used=new bool[COMM_OPT_MAX*block_num];
    memset(used,0,sizeof(bool)*COMM_OPT_MAX*block_num);
    ll bar[10];
    ll bar_t=1;
    for(int i=0;i<8;i++)
    {
        bar[i]=bar_t;
        bar_t<<=logmax;
    }
    for(int i=0;i<n;i++)
    {
            if(f[i]==0)
                continue;
            
            if(f[i]<0)
            {
                ll tmp=-f[i];
                for(int j=0;j<block_num;j++)
                {
                    if(tmp<bar[j])
                        break;
                    ll fnow=(tmp>>(logmax*j))&65535;
                    W[fnow+(j<<logmax)]-=g[i];
                    used[fnow+(j<<logmax)]=1;
                }
            }
            else
            {
                ll tmp=f[i];
                for(int j=0;j<block_num;j++)
                {
                    if(tmp<bar[j])
                        break;
                    ll fnow=(tmp>>(logmax*j))&65535;
                    W[fnow+(j<<logmax)]+=g[i];
                    used[fnow+(j<<logmax)]=1;
                }
            }
    }
    G1 gg[logmax*block_num];
    for(int j=0;j<logmax*block_num;j++)
        gg[j].clear();
    for(int j=0;j<COMM_OPT_MAX*block_num;j++)
    {
        if(used[j])
        {
            int jj=j%COMM_OPT_MAX;
            int blk=j/COMM_OPT_MAX;
            for(int k=0;k<logmax;k++)
            {
                if(jj&(1<<k))
                    gg[k+logmax*blk]+=W[j];
            }
            W[j].clear();
            used[j]=0;            
        }
    }
    for(int j=0;j<logmax*block_num;j++)
    {
        if(j>60)
        {
            G1 gd=gg[j]*(1ll<<48);
            ret+=gd*(1ll<<(j-48)); //split
        }
        else
            ret+=gg[j]*(1ll<<j);
    }
    delete []used;
    return ret;
}


G1 perdersen_commit(G1* g,int* f,int n,G1* W)
{
    G1 ret;
    ret.clear();
    bool *used=new bool[COMM_OPT_MAX];
    memset(used,0,sizeof(bool)*COMM_OPT_MAX);
    for(int i=0;i<n;i++)
    {
            if(f[i]==0)
                continue;
            
            if(f[i]<0)
            {
                W[-f[i]]-=g[i];
                used[-f[i]]=1;
                assert(-f[i]<COMM_OPT_MAX);
            }
            else
            {
                W[f[i]]+=g[i];
                used[f[i]]=1;
                assert(f[i]<COMM_OPT_MAX);
            }
    }
    //t.stop("add ",false);
    const int logn=log2(COMM_OPT_MAX)+1;
    G1 gg[40];
    for(int j=0;j<logn;j++)
        gg[j].clear();
    for(int j=1;j<COMM_OPT_MAX;j++)
    {
        if(used[j])
        {
            for(int k=0;k<logn;k++)
            {
                if(j&(1<<k))
                    gg[k]+=W[j];
            }
            W[j].clear();
            used[j]=0;            
        }
    }
    for(int j=0;j<logn;j++)
        ret+=gg[j]*(1<<j);

    //t.stop("accu",false);
    //t.stop("ALL: ",true);

    delete []used;
    return ret;
}

G1 perdersen_commit(G1* g,Fr* f,int n)
{
    G1 ret;
    G1::mulVec(ret,g,f,n);
    return ret;
}

Fr lagrange(Fr *r,int l,int k)
{
    assert(k>=0 && k<(1<<l));
    Fr ret=1;
    for(int i=0;i<l;i++)
    {
        if(k&(1<<i))
            ret*=r[i];
        else
            ret*=1-r[i];
    }
    return ret;
}
void brute_force_compute_LR(Fr* L,Fr* R,Fr* r,int l)
{
    int halfl=l/2,c=l-halfl;
    for(int k=0;k<(1<<c);k++)
        L[k]=lagrange(r,c,k);
    for(int k=0;k<(1<<halfl);k++)
        R[k]=lagrange(r+c,halfl,k);
}

Fr brute_force_compute_eval(Fr* w,Fr* r,int l)
{
    Fr ret=0;
    for(int k=0;k<(1<<l);k++)
        ret+=lagrange(r,l,k)*w[k];
    return ret;
}


G1 compute_Tprime(int l,Fr* R,G1* Tk) 
{
    //w has 2^l length
    //assert(l%2==0);
    int halfl=l/2;
    int rownum=(1<<halfl),colnum=(1<<(l-halfl));
    G1 ret=perdersen_commit(Tk,R,rownum);
    return ret;
}

G1 compute_RT_singlethread(Fr*w ,Fr*R,int l,G1*g,Fr*& ret) // L is row number length
{
    int halfl=l/2;
    int rownum=(1<<halfl),colnum=(1<<(l-halfl));
    Fr* res=new Fr[colnum];
    for(int i=0;i<colnum;i++)
        res[i]=0;
    timer TT;
    TT.start();
    for(int j=0;j<colnum;j++)
    for(int i=0;i<rownum;i++)
    {
        if(!w[j+i*colnum].isZero())
            res[j]+=R[i]*w[j+i*colnum];  // mat mult  (1,row)*(row,col)=(1,col)
    }
    TT.stop();
    G1 comm=perdersen_commit(g,res,colnum);
    ret=res;
    return comm;
}
G1 compute_RT(Fr *w, Fr *R, int l, G1 *g, Fr *&ret) // L is row number length
{
    int halfl = l / 2;
    int rownum = (1 << halfl), colnum = (1 << (l - halfl));
    Fr *res = new Fr[colnum];
    for (int i = 0; i < colnum; i++)
        res[i] = 0;

    auto start = std::chrono::high_resolution_clock::now();

    // Define the number of threads
    const int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    std::mutex res_mutex;

    auto worker = [&](int start_col, int end_col) {
        Fr *local_res = new Fr[colnum]; // Thread-local result storage
        for (int i = start_col; i < end_col; i++)
            local_res[i] = 0;

        for (int j = start_col; j < end_col; j++) {
            for (int i = 0; i < rownum; i++) {
                if (!w[j + i * colnum].isZero()) {
                    local_res[j] += R[i] * w[j + i * colnum];
                }
            }
        }

        // Merge local results into global result
        std::lock_guard<std::mutex> guard(res_mutex);
        for (int j = start_col; j < end_col; j++) {
            res[j] += local_res[j];
        }

        delete[] local_res; // Free thread-local storage
    };

    // Split columns into chunks for each thread
    int chunk_size = (colnum + num_threads - 1) / num_threads;
    for (int t = 0; t < num_threads; t++) {
        int start_col = t * chunk_size;
        int end_col = std::min(start_col + chunk_size, colnum);
        if (start_col < colnum) { // Avoid starting threads if no work is left
            threads.emplace_back(worker, start_col, end_col);
        }
    }

    // Join threads
    for (auto &th : threads) {
        if (th.joinable())
            th.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    // Commit the result
    G1 comm = perdersen_commit(g, res, colnum);
    ret = res;
    return comm;
}

G1 gen_gi(G1* g,int n)
{
    G1 base;
    for(int i=0;i<n;i++)
    {
        Fr tmp;
        tmp.setByCSPRNG();
        g[i]=base*tmp;
    }
    return base;
}

double blt_vtime=0;
Pack bullet_reduce(G1 gamma, Fr*a,G1*g,int n,G1& G,Fr* x,Fr y,bool need_free) // length n
{
    timer vtimer;
    if(n==1)
    {
        Pack p(gamma,a[0],g[0],x[0],y);
        return p;
    }
    
    //step2  prover fold
    G1 gamma_minus1,gamma_1;
    Fr x1a2=0,x2a1=0;
    for(int i=0;i<n/2;i++)
    {
        x1a2+=x[i]*a[n/2+i];
        x2a1+=x[n/2+i]*a[i];
    }
    gamma_minus1=G*x1a2+perdersen_commit(g+n/2,x,n/2);
    gamma_1=G*x2a1+perdersen_commit(g,x+n/2,n/2);
    Fr c,invc;
    c.setByCSPRNG();  // step3 V choose random c
    //prover verifier both comp
    vtimer.start();
    Fr::inv(invc,c);
    G1 gamma_prime=gamma_minus1*c*c+gamma_1*invc*invc+gamma;
    Fr* aprime=new Fr[n/2];       
    for(int i=0;i<n/2;i++)
        aprime[i]=a[i]*invc+a[i+n/2]*c;
    G1* gprime=new G1[n/2];        
    if(n<2048)
    {   
        for(int i=0;i<n/2;i++)
            gprime[i]=g[i]*invc+g[i+n/2]*c;
    }
    else  
    {
        int num_threads=16;
        int chunk_size = n/2/num_threads;
        std::vector<std::thread> threads1;
        auto worker1 = [&](int start_id, int end_id) 
        {
            for(int i=start_id;i<end_id;i++)
                gprime[i]=g[i]*invc+g[i+n/2]*c;
        };
        for (int t = 0; t < num_threads; t++) 
        {
            int start_col = t * chunk_size;
            int end_col = min(start_col + chunk_size, n/2);
            threads1.emplace_back(worker1, start_col, end_col);
        }
        for (auto &th : threads1) 
        {
            if (th.joinable())
                th.join();
        }
    }
    vtimer.stop();

    blt_vtime+=vtimer.elapse_sec();

    //prover single compute
    Fr* xprime=new Fr[n/2];         
    Fr yprime;
    for(int i=0;i<n/2;i++)
        xprime[i]=c*x[i]+invc*x[i+n/2];
    yprime=c*c*x1a2+invc*invc*x2a1+y;

    if(need_free)
    {
        delete []a;
        delete []g;
        delete []x;
    }
    
    return bullet_reduce(gamma_prime,aprime,gprime,n/2,G,xprime,yprime,true);
}   

void prove_dot_product(G1 comm_x, G1 comm_y, Fr* a, G1*g ,G1& G,Fr* x,Fr y,int n)  // y= <a,x> , 
{
    G1 gamma=comm_x+comm_y;
    Pack p=bullet_reduce(gamma,a,g,n,G,x,y);
    assert(p.y==p.x*p.a);
    assert(p.gamma==p.g*p.x+G*p.y);
}
static ThreadSafeQueue<int> workerq,endq;


void ll_commit_worker(G1*& Tk,G1*& g, ll*& row,int colnum,G1*& W)
{
    int idx;
    while (true)
    {
            bool ret=workerq.TryPop(idx);
            if(ret==false)
                return;
            Tk[idx]=perdersen_commit(g,row+idx*colnum,colnum,W);
            endq.Push(idx);
    }
}
G1* prover_commit(ll* w, G1* g, int l,int thread_n) //compute Tk, int version with pippenger
{
    int halfl=l/2;
    int rownum=(1<<halfl),colnum=(1<<(l-halfl));
    G1 *Tk=new G1[rownum];
    ll* row=new ll[1<<l];
    G1** W=new G1*[thread_n];
    for(int i=0;i<thread_n;i++)
        W[i]=new G1[COMM_OPT_MAX*block_num];
    for(int i=0;i<thread_n;i++)
        memset(W[i],0,sizeof(G1)*COMM_OPT_MAX*block_num);
    for (u64 i = 0; i < rownum; ++i)  //work for rownum 
        workerq.Push(i);
    for(int i=0;i<thread_n;i++)
    {
        thread t(ll_commit_worker,std::ref(Tk),std::ref(g),std::ref(w),colnum,std::ref(W[i])); 
        t.detach();
    }
    while(!workerq.Empty())
        this_thread::sleep_for (std::chrono::microseconds(10));
    while(endq.Size()!=rownum)
        this_thread::sleep_for (std::chrono::microseconds(10));
    endq.Clear();
    assert(endq.Size()==0);
    for(int i=0;i<thread_n;i++)
        delete [] W[i];
    delete []W;
    return Tk;
}



Fr prover_evaluate(Fr*w ,Fr*r,G1& G,G1* g, Fr*L,Fr*R,int l)  // nlogn brute force 
{
    int halfl=l/2;
    int rownum=(1<<halfl),colnum=(1<<(l-halfl));
    timer t(true);
    t.start();
    brute_force_compute_LR(L,R,r,l);
    Fr eval=brute_force_compute_eval(w,r,l);
    t.stop("eval total ",true,false);
    return eval;
}
namespace hyrax
{
pair<double,double> open(Fr*w,Fr*r,Fr eval,G1&G,G1*g,Fr*L,Fr*R,G1*tk,int l)
{
    double prover_time=0,verifier_time=0;
    int halfl=l/2;
    int rownum=(1<<halfl),colnum=(1<<(l-halfl));
    timer verf;
    verf.start();
    Fr* RT=new Fr[colnum];
    compute_RT(w,R,l,g,RT);  
    verf.stop();
    prover_time+=verf.elapse_sec();
    verf.start();
    G1 tprime=compute_Tprime(l,R,tk);
    verf.stop();
    prover_time+=verf.elapse_sec();
    verifier_time+=verf.elapse_sec();
    verf.start();
    prove_dot_product(tprime, G*eval, L, g , G,RT,eval,colnum);
    verf.stop();
    prover_time+=verf.elapse_sec();
    verifier_time+=blt_vtime;
    return make_pair(prover_time,verifier_time);
}
}