#ifndef HYRAX_RP_DEFINE
#define HYRAX_RP_DEFINE
#include <iostream>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>         // std::thread
#include <mcl/bn256.hpp>
#include "timer.hpp"
#include "typedef.hpp"
#include "hyrax_util.hpp"
using namespace std;
using namespace mcl::bn;

G1 range_proof_perdersen_commit(G1* g,ll* f,int n,G1* W=NULL); //support 2^80, optimized using pippenger
G1 range_proof_perdersen_commit_classic(G1* g,Fr* f,int n);
G1 range_proof_perdersen_commit_redundant(G1* g,Fr* f,ll* id,int n,int m); 
Fr range_proof_lagrange(Fr *r,int l,int k);
G1* range_proof_prover_commit_fr(ll* w, Fr* f,int m,G1* g, int l,int thread_n);
G1* range_proof_prover_commit_fr_general(Fr* w, G1* g, int l,int thread_n);
Fr* range_proof_get_eq(Fr*r, int l);

G1 range_proof_gen_gi(G1* g,int n);

bool range_proof_prove_dot_product(G1 comm_x, G1 comm_y, Fr* a, Fr* x,Fr y, G1*g ,G1& G,int n);

G1* range_proof_prover_commit(ll* w, G1* g, int l,int thread_n=1);
Fr range_proof_prover_evaluate(ll*w ,Fr*r,int l);  
void range_proof_open(ll*w,Fr*r,Fr eval,G1&G,G1*g,G1*comm,int l);
void range_proof_open(Fr*w,Fr*r,Fr eval,G1&G,G1*g,G1*comm,int l);

#endif
