#ifndef ZKCNN_RANGE_PROVER_HPP
#define ZKCNN_RANGE_PROVER_HPP

#include "global_var.hpp"
#include "circuit.h"
#include "polynomial.h"
#include "hyrax_rp.hpp"
using std::unique_ptr;
const int MAXL=24;  

class neuralNetwork;
enum class NonlinearOpType {
    Attention, LayerNorm, GELU, NonLinear
};

inline int next_power_of_two_exp(int n) {
    if (n <= 1) return 0;
    return 32 - __builtin_clz(n - 1);
}

class Constraint {
public:
    int query_size;
    int range_size;
    ll * inputs;
    void prepare();
};

class OP {
public:
    NonlinearOpType op_type;
    vector<Constraint> constraints;
};

class range_prover {
public:
    // Constructor
    range_prover(int layer_num = 0, int head_num = 0, int head_dim = 0, 
                int attn_dim = 0, int linear_dim = 0, int seq_len = 0, 
                int threads = 32, int merged = 0) 
        : LayerNum(layer_num), HeadNum(head_num), HeadDim(head_dim),
          AttnDim(attn_dim), LinearDim(linear_dim), seq_len(seq_len), 
          thread_num(threads), merged(merged) {}

    struct SC_Return {
        Fr* random;
        Fr claim_f;
        Fr claim_g;
    };
    int merged = 0;
    SC_Return sumcheck_deg1(int l, Fr* f, Fr S);
    SC_Return sumcheck_deg3(int l, Fr* r, Fr* f, Fr* g, Fr S);
    int LayerNum, HeadNum, HeadDim, AttnDim, LinearDim, seq_len;// llm related parameters
    int thread_num = 32;
    vector<OP> ops;
    G1 g[1<<(MAXL/2)],GG;
    Fr r[MAXL];
    ll * inputs = new ll[1<<MAXL];
    void init();
    void push_back(NonlinearOpType op_type, const std::vector<std::pair<int, int>>& constraint_params);
    timer prove_timer;
    timer prepare_timer;
    void range_prove(ll * inputs,int range,int query_size,int thread_num);
    void logup(ll * f,ll * t,int m,int n,int thread_num);
    double prove();
    void build();
private:
    
};
__int128 convert(Fr x)	;

#endif //ZKCNN_RANGE_PROVER_HPP
