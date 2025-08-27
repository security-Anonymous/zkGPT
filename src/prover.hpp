#ifndef ZKCNN_PROVER_HPP
#define ZKCNN_PROVER_HPP

#include "global_var.hpp"
#include "circuit.h"
#include "polynomial.h"
#include "hyrax.hpp"
using std::unique_ptr;

class neuralNetwork;
struct Commit_return
{
    int l;
    Fr* w;
    ll* ww;
    G1 G,*g;
    G1* comm;
};

class prover {
public:
    int fc_row[100],fc_col[100],fc_start_id[100], fc_input_row[100],fc_input_col[100],fc_input_id[100];
    int fc_real_row[100],fc_real_col[100],fc_real_input_row[100],fc_real_input_col[100]; // padded
    u32 total[2], total_size[2];
    vector<linear_poly> mult_array[2];
    vector<linear_poly> V_mult[2];
    vector<linear_poly> tmp_mult_array[2];
    vector<linear_poly> tmp_V_mult[2];
    int** mat_val;
    void init();

    void sumcheckInitAll(const vector<F>::const_iterator &r_0_from_v);
    void sumcheckInit(const F &alpha_0, const F &beta_0);
    void sumcheckDotProdInitPhase1();
    void sumcheckInitPhase1(const F &relu_rou_0);
    void sumcheckInitPhase2();

    cubic_poly sumcheckDotProdUpdate1(const F &previous_random);
    quadratic_poly sumcheckUpdate1(const F &previous_random);
    quadratic_poly sumcheckUpdate2(const F &previous_random);

    F Vres(const vector<F>::const_iterator &r, u32 output_size, u8 r_size,int layer_id,int start=0) ;
    vector<F> lasso_mult_v,lasso_v_mult;
    void sumcheckDotProdFinalize1(const F &previous_random, F &claim_1);
    void sumcheckFinalize1(const F &previous_random, F &claim_0, F &claim_1);
    void sumcheckFinalize2(const F &previous_random, F &claim_0, F &claim_1);
    void sumcheck_lasso_Finalize(const F &previous_random, F &claim_1);

    void sumcheckLassoInit(const vector<F> &s_u, const vector<F> &s_v,const vector<vector<F>>& r_uu, const vector<vector<F>>& r_vv);
    quadratic_poly sumcheckLassoUpdate(const F &previous_random);
    quadratic_poly sumcheckUpdateEach_Lasso(const F &previous_random, bool idx) ;
    void commitInput(const vector<G1> &gens,int thread_n=1);

    timer prove_timer;
    vector<double> throw_time[500];
    double proveTime() const{ return 0;}// { return prove_timer.elapse_sec(); }
    double proofSize() const { return 0;}// (double) proof_size / 1024.0; }
    double polyProverTime() const { return 0;}//poly_p -> getPT(); }
    double polyProofSize() const { return 0;}//{ return poly_p -> getPS(); }
    
    layeredCircuit C;
    vector<vector<F>> val;        // the output of each gate
    F getCirValue(u8 layer_id, const vector<u32> &ori, u32 u);
    u64 proof_size;
    vector<G> gens;
    Commit_return cc;
    
private:
    quadratic_poly sumcheckUpdateEach(const F &previous_random, bool idx);
    quadratic_poly sumcheckUpdate(const F &previous_random, vector<F> &r_arr);
    

    vector<F>::iterator r_0, r_1;         // current positions
    vector<vector<F>> r_u, r_v;             // next positions

    vector<F> beta_g;
    
    F add_term;
    

    F V_u0, V_u1;

    F alpha, beta, relu_rou;

    

    
    u8 round;          // step within a sumcheck
    u8 sumcheck_id;    // the level

    friend neuralNetwork;
};
__int128 convert(Fr x)	;

#endif //ZKCNN_PROVER_HPP
