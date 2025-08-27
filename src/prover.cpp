#include "prover.hpp"
#include <iostream>
#include <utils.hpp>

static vector<F> beta_gs, beta_u;
using namespace mcl::bn;
using std::unique_ptr;

linear_poly interpolate(const F &zero_v, const F &one_v) 
{
    return {one_v - zero_v, zero_v};
}

F prover::getCirValue(u8 layer_id, const vector<u32> &ori, u32 u) {
    return !layer_id ? val[0][ori[u]] : val[layer_id][u];
}

void prover::init() 
{
    proof_size = 0;
    r_u.resize(C.size + 1);
    r_v.resize(C.size + 1);
    const int SIZE=28;
    V_mult[0].resize(1<<SIZE);
    V_mult[1].resize(1<<SIZE);
    mult_array[0].resize(1<<SIZE);
    mult_array[1].resize(1<<SIZE);
    tmp_V_mult[0].resize(1<<SIZE);
    tmp_V_mult[1].resize(1<<SIZE);
    tmp_mult_array[0].resize(1<<SIZE);
    tmp_mult_array[1].resize(1<<SIZE);
    lasso_mult_v.resize(1<<SIZE);
    for(int i=0;i<(1<<SIZE);i++)
        lasso_mult_v[i]=0;
}

/**
 * This is to initialize all process.
 *
 * @param the random point to be evaluated at the output layer
 */
void prover::sumcheckInitAll(const vector<F>::const_iterator &r_0_from_v) 
{
    sumcheck_id = C.size;
    i8 last_bl = C.circuit[sumcheck_id - 1].bit_length;
    r_u[sumcheck_id].resize(last_bl);
    prove_timer.start();
    for (int i = 0; i < last_bl; ++i) 
        r_u[sumcheck_id][i] = r_0_from_v[i];
    prove_timer.stop();
}

/**
 * This is to initialize before the process of a single layer.
 *
 * @param the random combination coefficiants for multiple reduction points
 */
void prover::sumcheckInit(const F &alpha_0, const F &beta_0) 
{
    prove_timer.start();
    auto &cur = C.circuit[sumcheck_id];
    alpha = alpha_0;
    beta = beta_0;
    r_0 = r_u[sumcheck_id].begin();
    r_1 = r_v[sumcheck_id].begin();
    --sumcheck_id;
    prove_timer.stop();
}
static ThreadSafeQueue<int> workerq,endq;

void sc_phase1_uni_worker(vector<uniGate> &beg, std::vector<linear_poly> (&mult_array)[2],vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];

        for (size_t i = l; i <r ; ++i) 
        {
            auto &gate = beg[i];
            bool idx = gate.lu != 0;
            mult_array[idx][gate.u] = mult_array[idx][gate.u] + beta_g[gate.g] * gate.sc;
        }
        endq.Push(idx);
    }
}

void sc_phase1_bin_worker(layer& cur, vector<binGate> &beg, std::vector<linear_poly> (&mult_array)[2],F& V_u0,F&V_u1,vector<vector<F> >& val,vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R,u8 sumcheck_id) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];

        for (size_t i = l; i <r ; ++i) 
        {
            auto &gate = beg[i];
            bool idx = gate.getLayerIdU(sumcheck_id) != 0;
            auto val_lv =  !gate.getLayerIdV(sumcheck_id) ? val[0][cur.ori_id_v[gate.v]] : val[gate.getLayerIdV(sumcheck_id)][gate.v];
            mult_array[idx][gate.u] = mult_array[idx][gate.u] + val_lv * beta_g[gate.g] * gate.sc;  // Ahg for phase 1
        }
        endq.Push(idx);
    }
}

void prover::sumcheckInitPhase1(const F &relu_rou_0) 
{
    //fprintf(stderr, "sumcheck level %d, phase1 init start\n", sumcheck_id);
    auto &cur = C.circuit[sumcheck_id];
    total[0] = ~cur.bit_length_u[0] ? 1ULL << cur.bit_length_u[0] : 0;
    total_size[0] = cur.size_u[0];
    total[1] = ~cur.bit_length_u[1] ? 1ULL << cur.bit_length_u[1] : 0;
    total_size[1] = cur.size_u[1];

    r_u[sumcheck_id].resize(cur.max_bl_u);
    timer useless_t;
    useless_t.start();
    beta_g.resize(1ULL << cur.bit_length);
    relu_rou = relu_rou_0;
    add_term.clear();
    
    for (int b = 0; b < 2; ++b)
        for (u32 u = 0; u < total[b]; ++u)
            mult_array[b][u].clear();
    
    
        for (int b = 0; b < 2; ++b) 
        {
            auto dep = !b ? 0 : sumcheck_id - 1;
            for (u32 u = 0; u < total[b]; ++u) 
            {
                if (u >= cur.size_u[b])
                    V_mult[b][u].clear();
                else 
                {
                    V_mult[b][u] = getCirValue(dep, cur.ori_id_u, u);
                }
            }
        }
    useless_t.stop();
    throw_time[sumcheck_id].push_back(useless_t.elapse_sec());
    prove_timer.start();
    
    initBetaTable(beta_g, cur.bit_length, r_0, r_1, alpha, beta);
    
    if(cur.uni_interval.size()>=2)
    {
        const int thd=32;
        int *L=new int [cur.uni_interval.size()],*R=new int [cur.uni_interval.size()];
        for (u64 j = 0; j <cur.uni_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur.uni_interval[j].first;
                R[j]=cur.uni_interval[j].second;
        }
        for(int i=0;i<thd;i++)
        {
            thread t(sc_phase1_uni_worker, std::ref(cur.uni_gates),std::ref(mult_array),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R)); 
            t.detach();
        }
        while(!workerq.Empty())
            this_thread::sleep_for (std::chrono::microseconds(1));
        while(endq.Size()!=cur.uni_interval.size())
            this_thread::sleep_for (std::chrono::microseconds(1));
        endq.Clear();
    }
    else for (auto &gate: cur.uni_gates) 
        {
            bool idx = gate.lu != 0;
            mult_array[idx][gate.u] = mult_array[idx][gate.u] + beta_g[gate.g] * gate.sc;
        }
    if(cur.bin_interval.size()>=2)
    {
        const int thd=32;
        int *L=new int [cur.bin_interval.size()],*R=new int [cur.bin_interval.size()];
        for (u64 j = 0; j <cur.bin_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur.bin_interval[j].first;
                R[j]=cur.bin_interval[j].second;
        }
        for(int i=0;i<thd;i++)
        {
            thread t(sc_phase1_bin_worker, std::ref(cur),std::ref(cur.bin_gates),std::ref(mult_array),std::ref(V_u0),std::ref(V_u1),std::ref(val),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R),sumcheck_id); 
            t.detach();
        }
        while(!workerq.Empty())
            this_thread::sleep_for (std::chrono::microseconds(1));
        while(endq.Size()!=cur.bin_interval.size())
            this_thread::sleep_for (std::chrono::microseconds(1));
        endq.Clear();
    }
    else  for (auto &gate: cur.bin_gates) 
        {
            bool idx = gate.getLayerIdU(sumcheck_id) != 0;
            auto val_lv = getCirValue(gate.getLayerIdV(sumcheck_id), cur.ori_id_v, gate.v);
            mult_array[idx][gate.u] = mult_array[idx][gate.u] + val_lv * beta_g[gate.g] * gate.sc;  // Ahg for phase 1
        }
    round = 0;
    prove_timer.stop();
    //fprintf(stderr, "sumcheck level %d, phase1 init finished\n", sumcheck_id);
}


void sc_phase2_uni_worker( vector<uniGate> &beg, F& sum_value,F& V_u0,F&V_u1,vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        F ss;
        ss.clear();
        for (size_t i = l; i <r ; ++i) 
        {
            auto &gate = beg[i];
            auto V_u = !gate.lu ? V_u0 : V_u1;                  //V_u0 is claim 0, V_u1 is claim 1
            ss +=beta_g[gate.g] * beta_u[gate.u] * V_u * gate.sc;
        }
        sum_value+=ss;
        endq.Push(idx);
    }
}

void sc_phase2_bin_worker( vector<binGate> &beg, std::vector<linear_poly> (&mult_array)[2],F& V_u0,F&V_u1,vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R,u8 sumcheck_id) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];

        for (size_t i = l; i <r ; ++i) 
        {
            auto &gate = beg[i];
            bool idx = gate.getLayerIdV(sumcheck_id);
            auto V_u = !gate.getLayerIdU(sumcheck_id) ? V_u0 : V_u1;
            mult_array[idx][gate.v] =mult_array[idx][gate.v]+beta_g[gate.g] * beta_u[gate.u] * V_u * gate.sc;
        }
        endq.Push(idx);
    }
}
void prover::sumcheckInitPhase2() 
{
    //fprintf(stderr, "sumcheck level %d, phase2 init start\n", sumcheck_id);
    auto &cur = C.circuit[sumcheck_id];
    total[0] = ~cur.bit_length_v[0] ? 1ULL << cur.bit_length_v[0] : 0;
    total_size[0] = cur.size_v[0];
    total[1] = ~cur.bit_length_v[1] ? 1ULL << cur.bit_length_v[1] : 0;
    total_size[1] = cur.size_v[1];
    i8 fft_bl = cur.fft_bit_length;
    i8 cnt_bl = cur.max_bl_v;

    timer useless_time;
    useless_time.start();
    r_v[sumcheck_id].resize(cur.max_bl_v);


    beta_u.resize(1ULL << cur.max_bl_u);

    

    add_term.clear();
    for (int b = 0; b < 2; ++b) 
    {
        for (u32 v = 0; v < total[b]; ++v)
            mult_array[b][v].clear();
    }
    useless_time.stop();
    throw_time[sumcheck_id].push_back(useless_time.elapse_sec());
    prove_timer.start();
    initBetaTable(beta_u, cur.max_bl_u, r_u[sumcheck_id].begin(), F_ONE,32); //  beta_u is U in the code
    for (int b = 0; b < 2; ++b) 
    {
        auto dep = !b ? 0 : sumcheck_id - 1;
        for (u32 v = 0; v < total[b]; ++v) 
        {
            V_mult[b][v] = v >= cur.size_v[b] ? F_ZERO : getCirValue(dep, cur.ori_id_v, v);
        }
    }
    
    if(cur.uni_interval.size()>=2)
    {
        const int thd=32;
        F sum[40];
        
        int *L=new int [cur.uni_interval.size()],*R=new int [cur.uni_interval.size()];
        for (u64 j = 0; j <cur.uni_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur.uni_interval[j].first;
                R[j]=cur.uni_interval[j].second;
        }
            for(int i=0;i<thd;i++)
            {
                sum[i].clear(); 
                thread t(sc_phase2_uni_worker, std::ref(cur.uni_gates),std::ref(sum[i]),std::ref(V_u0),std::ref(V_u1),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R)); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(1));
            while(endq.Size()!=cur.uni_interval.size())
                this_thread::sleep_for (std::chrono::microseconds(1));
            endq.Clear();
            for(int i=0;i<thd;i++)
                add_term+=sum[i];
    }
    else for (auto &gate: cur.uni_gates) 
    {
        auto V_u = !gate.lu ? V_u0 : V_u1;                  //V_u0 is claim 0, V_u1 is claim 1
        add_term = add_term + beta_g[gate.g] * beta_u[gate.u] * V_u * gate.sc;
    }
    if(cur.bin_interval.size()>=2)
    {
        const int thd=32;
        int *L=new int [cur.bin_interval.size()],*R=new int [cur.bin_interval.size()];
        for (u64 j = 0; j <cur.bin_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur.bin_interval[j].first;
                R[j]=cur.bin_interval[j].second;
        }
            for(int i=0;i<thd;i++)
            {
                thread t(sc_phase2_bin_worker, std::ref(cur.bin_gates),std::ref(mult_array),std::ref(V_u0),std::ref(V_u1),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R),sumcheck_id); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(1));
            while(endq.Size()!=cur.bin_interval.size())
                this_thread::sleep_for (std::chrono::microseconds(1));
            endq.Clear();
    }
    else for (auto &gate: cur.bin_gates) 
    {
        bool idx = gate.getLayerIdV(sumcheck_id);
        auto V_u = !gate.getLayerIdU(sumcheck_id) ? V_u0 : V_u1;
        mult_array[idx][gate.v] = mult_array[idx][gate.v] + beta_g[gate.g] * beta_u[gate.u] * V_u * gate.sc;
    }

    round = 0;
    prove_timer.stop();
}

void prover::sumcheckLassoInit(const vector<F> &s_u, const vector<F> &s_v,const vector<vector<F>>& r_uu, const vector<vector<F>>& r_vv) 
{
    
    sumcheck_id = 0;
    total[1] = (1ULL << C.circuit[sumcheck_id].bit_length);
    total_size[1] = C.circuit[sumcheck_id].size;

    r_u[0].resize(C.circuit[0].bit_length);
    timer ggg;
    ggg.start();

    i8 max_bl = 0;
    for (int i = sumcheck_id + 1; i < C.size; ++i)
        max_bl = max(max_bl, max(C.circuit[i].bit_length_u[0], C.circuit[i].bit_length_v[0]));
    beta_g.resize(1ULL << max_bl);
    for (u8 i = sumcheck_id + 1; i < C.size; ++i) 
    {
        i8 bit_length_i = C.circuit[i].bit_length_u[0];
        u32 size_i = C.circuit[i].size_u[0];
        //timer a,b;
        if (~bit_length_i) 
        {
            r_u[i].resize(C.circuit[i].max_bl_u);
            for(int j=0;j<C.circuit[i].max_bl_u;j++)
                r_u[i][j]=r_uu[i][j];
            initBetaTable(beta_g, bit_length_i, r_u[i].begin(), s_u[i - 1],32);
            for (u32 hu = 0; hu < size_i; ++hu) 
            {
                u32 u = C.circuit[i].ori_id_u[hu];
                lasso_mult_v[u] += beta_g[hu];
            }
        }
        bit_length_i = C.circuit[i].bit_length_v[0];
        size_i = C.circuit[i].size_v[0];
        if (~bit_length_i) 
        {
            r_v[i].resize(C.circuit[i].max_bl_v);
            for(int j=0;j<C.circuit[i].max_bl_v;j++)
                r_v[i][j]=r_vv[i][j];
            initBetaTable( beta_g, bit_length_i, r_v[i].begin(), s_v[i - 1],32);
            for (u32 hv = 0; hv < size_i; ++hv) 
            {
                u32 v = C.circuit[i].ori_id_v[hv];
                lasso_mult_v[v] += beta_g[hv];
            }
        }
    }
    round = 0;
    prove_timer.stop();
}

quadratic_poly prover::sumcheckUpdate1(const F &previous_random) {
    return sumcheckUpdate(previous_random, r_u[sumcheck_id]);
}

quadratic_poly prover::sumcheckUpdate2(const F &previous_random) {
    return sumcheckUpdate(previous_random, r_v[sumcheck_id]);
}

quadratic_poly prover::sumcheckUpdate(const F &previous_random, vector<F> &r_arr) 
{
    prove_timer.start();

    if (round) r_arr.at(round - 1) = previous_random;
    ++round;
    quadratic_poly ret;

    add_term = add_term * (F_ONE - previous_random);
    for (int b = 0; b < 2; ++b)
        ret = ret + sumcheckUpdateEach(previous_random, b);
    ret = ret + quadratic_poly(F_ZERO, -add_term, add_term);

    prove_timer.stop();
    proof_size += F_BYTE_SIZE * 3;
    return ret;
}



void sumcheckUpdate_worker(quadratic_poly &sum,vector<linear_poly> &tmp_v,vector<linear_poly> &tmp_mult,vector<linear_poly> &tmp_v_2,vector<linear_poly> &tmp_mult_2,int*&L,int*&R,Fr previous_random,int total_size)
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        for(int i=l;i<r;i++)
        {
            u32 g0 = i << 1, g1 = i << 1 | 1;
            if (g0 >= total_size) 
                break;
            tmp_v_2[i] = interpolate(tmp_v[g0].eval(previous_random), tmp_v[g1].eval(previous_random));
            tmp_mult_2[i] = interpolate(tmp_mult[g0].eval(previous_random), tmp_mult[g1].eval(previous_random));
            sum = sum + tmp_mult_2[i] * tmp_v_2[i];
        }
        endq.Push(idx);
    }
}
quadratic_poly prover::sumcheckUpdateEach(const F &previous_random, bool idx) 
{
    auto &tmp_mult = mult_array[idx];
    auto &tmp_v = V_mult[idx];
    auto &tmp_mult_2 = tmp_mult_array[idx];
    auto &tmp_v_2 = tmp_V_mult[idx];

    if (total[idx] == 1) 
    {
        tmp_v[0] = tmp_v[0].eval(previous_random);
        tmp_mult[0] = tmp_mult[0].eval(previous_random);
        add_term = add_term + tmp_v[0].b * tmp_mult[0].b;
    }

    quadratic_poly ret;
    ret.clear();
    if(total[idx]<(1<<15))
    {
        for (u32 i = 0; i < (total[idx] >> 1); ++i) 
        {
            u32 g0 = i << 1, g1 = i << 1 | 1;
            if (g0 >= total_size[idx]) 
            {
                tmp_v[i].clear();
                tmp_mult[i].clear();
                continue;
            }
            tmp_v[i] = interpolate(tmp_v[g0].eval(previous_random), tmp_v[g1].eval(previous_random));
            tmp_mult[i] = interpolate(tmp_mult[g0].eval(previous_random), tmp_mult[g1].eval(previous_random));
            ret = ret + tmp_mult[i] * tmp_v[i];
        }
    }
    else
    {
        timer tt_f1,tt_f2;
        tt_f1.start();
        const int k=10;
        int total_work=(total[idx] >> 1);
        int *L=new int [(1<<k)],*R=new int [(1<<k)];
        const int thd=32;
        for (u64 j = 0; j < (1<<k); ++j) 
        {
            workerq.Push(j);
            L[j]=(total_work>>k)*j;
            R[j]=(total_work>>k)*(1+j);
        }
        quadratic_poly qp[thd];
        for(int j=0;j<thd;j++)
        {
            thread t(sumcheckUpdate_worker,std::ref(qp[j]),std::ref(tmp_v),std::ref(tmp_mult),std::ref(tmp_v_2),std::ref(tmp_mult_2),std::ref(L),std::ref(R),previous_random,total_size[idx]); 
            t.detach();
        }
        while(endq.Size()!=(1<<k))
            this_thread::sleep_for(std::chrono::microseconds(1));
        endq.Clear();
        for(int j=0;j<thd;j++)
            ret=ret+qp[j];
        tt_f1.stop();
        tt_f2.start();
        for(int i=0;i<total_work;i++)
        {
            if((i<<1)<total_size[idx])
            {
                tmp_mult[i]=tmp_mult_2[i];
                tmp_v[i]=tmp_v_2[i];
            }
            else
            {
                tmp_mult[i].clear();
                tmp_v[i].clear();
            }
        }
        tt_f2.stop();
    }
    
    total[idx] >>= 1;
    total_size[idx] = (total_size[idx] + 1) >> 1;

    return ret;
}


/**
 * This is to evaluate a multi-linear extension at a random point.
 *
 * @param the value of the array & random point & the size of the array & the size of the random point
 * @return sum of `values`, or 0.0 if `values` is empty.
 */
F prover::Vres(const vector<F>::const_iterator &r, u32 output_size, u8 r_size,int layer_id,int start) 
{
    prove_timer.start();

    vector<F> output(output_size);
    for (u32 i = 0; i < output_size; ++i)
        output[i] = val[layer_id][i+start];
    u32 whole = 1ULL << r_size;
    for (u8 i = 0; i < r_size; ++i) {
        for (u32 j = 0; j < (whole >> 1); ++j) {
            if (j > 0)
                output[j].clear();
            if ((j << 1) < output_size)
                output[j] = output[j << 1] * (F_ONE - r[i]);
            if ((j << 1 | 1) < output_size)
                output[j] = output[j] + output[j << 1 | 1] * (r[i]);
        }
        whole >>= 1;
    }
    F res = output[0];

    prove_timer.stop();
    proof_size += F_BYTE_SIZE;
    return res;
}

void prover::sumcheckFinalize1(const F &previous_random, F &claim_0, F &claim_1) {
    prove_timer.start();
    r_u[sumcheck_id].at(round - 1) = previous_random;
    V_u0 = claim_0 = total[0] ? V_mult[0][0].eval(previous_random) : (~C.circuit[sumcheck_id].bit_length_u[0]) ? V_mult[0][0].b : F_ZERO;
    V_u1 = claim_1 = total[1] ? V_mult[1][0].eval(previous_random) : (~C.circuit[sumcheck_id].bit_length_u[1]) ? V_mult[1][0].b : F_ZERO;
    prove_timer.stop();

    mult_array[0].clear();
    mult_array[1].clear();
    V_mult[0].clear();
    V_mult[1].clear();
    proof_size += F_BYTE_SIZE * 2;
}

void prover::sumcheckFinalize2(const F &previous_random, F &claim_0, F &claim_1) {
    prove_timer.start();
    r_v[sumcheck_id].at(round - 1) = previous_random;
    claim_0 = total[0] ? V_mult[0][0].eval(previous_random) : (~C.circuit[sumcheck_id].bit_length_v[0]) ? V_mult[0][0].b : F_ZERO;
    claim_1 = total[1] ? V_mult[1][0].eval(previous_random) : (~C.circuit[sumcheck_id].bit_length_v[1]) ? V_mult[1][0].b : F_ZERO;
    prove_timer.stop();

    mult_array[0].clear();
    mult_array[1].clear();
    V_mult[0].clear();
    V_mult[1].clear();
    proof_size += F_BYTE_SIZE * 2;
}

void prover::sumcheck_lasso_Finalize(const F &previous_random, F &claim_1) 
{
    prove_timer.start();
    r_u[sumcheck_id].at(round - 1) = previous_random;
    claim_1 = total[1] ? V_mult[1][0].eval(previous_random) : V_mult[1][0].b;
    prove_timer.stop();
    proof_size += F_BYTE_SIZE;
}



void prover::commitInput(const vector<G1> &gens,int thr) 
{
    int len;
    if (C.circuit[0].size != (1ULL << C.circuit[0].bit_length)) 
    {
        len=val[0].size();
        val[0].resize(1ULL << C.circuit[0].bit_length);
        for (int i = C.circuit[0].size; i < val[0].size(); ++i)
            val[0][i].clear();
    }
    
    int l=ceil(log2(val[0].size()));
    ll* vi=new ll[1<<l];
    memset(vi,0,sizeof(ll)*(1<<l));
    ll mx=-1e9,mn=1e9;
    for(int i=0;i<val[0].size();i++)
    {
        vi[i]=convert(val[0][i]);
        mx=max(mx,vi[i]);
        mn=min(mn,vi[i]);
        
    }

    
    
    Fr* dat=new Fr[1<<l];
    memset(dat,0,sizeof(Fr)*(1<<l));
    memcpy(dat,val[0].data(),sizeof(Fr)*val[0].size());
    G1* ret=prover_commit(vi,(G1*)gens.data(),l,thr);
    cc.comm=ret;
    cc.G=gens.back();
    cc.g=(G1*)gens.data();
    cc.l=l;
    cc.w=dat;
    cc.ww=vi;
    val[0].resize(len);
}

__int128 convert(Fr x)	
{	
    int sign=0;	
    Fr abs;	
    if(x.isNegative())	
    {	
        sign=1;	
        abs=-x;	
    }	
    else	
        abs=x;	

    uint8_t bf[16]={0};	 //64 bit
    int size=abs.getLittleEndian(bf,16);	
    ll V=0;	
    for(int j=size-1;j>=0;j--)	
        V=V*256+bf[j];	
    if(sign)	
        V=-V;	
    return V;	
}