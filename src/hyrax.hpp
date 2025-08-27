#ifndef HYRAX_DEFINE
#define HYRAX_DEFINE
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>         // std::thread
#include <mcl/bn256.hpp>
#include "typedef.hpp"
#include "timer.hpp"
#include "hyrax_util.hpp"
using namespace std;
using namespace mcl::bn;
struct Hyrax_first_round
{
    //transmitted data
    Fr* tk;
    G1 comm_w;
};
struct Hyrax_second_round
{  
    //public data
    G1 *g, G;
    Fr* R;
    //prover only data
    Fr* LT;
    Fr eval;
};

G1 perdersen_commit(G1* g,ll* f,int n,G1* W=NULL); //support 2^80, optimized using pippenger
Fr lagrange(Fr *r,int l,int k);
void brute_force_compute_LR(Fr* L,Fr* R,Fr* r,int l);
Fr brute_force_compute_eval(Fr* w,Fr* r,int l);
G1 compute_Tprime(int l,Fr* L,G1* Tk) ;
G1 compute_LT(Fr*w ,Fr*L,int l,G1*g,Fr*& ret);
G1 gen_gi(G1* g,int n);
Pack bullet_reduce(G1 gamma, Fr*a,G1*g,int n,G1& G,Fr* x,Fr y,bool need_free=false);
void prove_dot_product(G1 comm_x, G1 comm_y, Fr* a, G1*g ,G1& G,Fr* x,Fr y,int n);
G1* prover_commit(ll* w, G1* g, int l,int thread_n=1);
G1* prover_commit(Fr* w, G1* g, int l,int thread=1);
G1* prover_commit(int* w, G1* g, int l,int thread=1);
Fr prover_evaluate(Fr*w ,Fr*r,G1& G,G1* g, Fr*L,Fr*R,int l);  // nlogn brute force 
namespace hyrax
{
pair<double,double> open(Fr*w,Fr*r,Fr eval,G1&G,G1*g,Fr*L,Fr*R,G1*tk,int l);
}


#endif
