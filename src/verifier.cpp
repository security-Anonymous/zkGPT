//
// Created by 69029 on 3/9/2021.
//
#include "verifier.hpp"
#include "global_var.hpp"
#include <utils.hpp>
#include <circuit.h>
#include <iostream>

using namespace std;
using namespace mcl::bn;
vector<F> beta_v;
static vector<F> beta_u, beta_gs;



verifier::verifier(prover *pr, const layeredCircuit &cir):
    p(pr), C(cir) 
{
    comm=p->cc;
    final_claim_u0.resize(C.size + 2);
    final_claim_v0.resize(C.size + 2);
    prover_time=verifier_time=0;
    matrix_time=0;
    r_u.resize(C.size + 2);
    r_v.resize(C.size + 2);
    // make the prover ready
    p->init();
}

F verifier::getFinalValue(const F &claim_u0, const F &claim_u1, const F &claim_v0, const F &claim_v1) {
    auto test_value = bin_value[0] * (claim_u0 * claim_v0)
                      + bin_value[1] * (claim_u1 * claim_v1)
                      + bin_value[2] * (claim_u1 * claim_v0)
                      + uni_value[0] * claim_u0
                      + uni_value[1] * claim_u1;

    return test_value;
}

void verifier::betaInitPhase1(u8 depth, const F &alpha, const F &beta, const vector<F>::const_iterator &r_0, const vector<F>::const_iterator &r_1, const F &relu_rou) {
    i8 bl = C.circuit[depth].bit_length;
    i8 fft_bl = C.circuit[depth].fft_bit_length;
    i8 fft_blh = C.circuit[depth].fft_bit_length - 1;
    i8 cnt_bl = bl - fft_bl, cnt_bl2 = C.circuit[depth].max_bl_u - fft_bl;


        
            beta_g.resize(1ULL << bl);
            initBetaTable(beta_g, C.circuit[depth].bit_length, r_0, r_1, alpha, beta );
            if (C.circuit[depth].zero_start_id < C.circuit[depth].size)
                for (u32 g = C.circuit[depth].zero_start_id; g < 1ULL << C.circuit[depth].bit_length; ++g)
                    beta_g[g] = beta_g[g] * relu_rou;
            beta_u.resize(1ULL << C.circuit[depth].max_bl_u);
            initBetaTable(beta_u, C.circuit[depth].max_bl_u, r_u[depth].begin(), F_ONE);
    
}

void verifier::betaInitPhase2(u8 depth) {
    beta_v.resize(1ULL << C.circuit[depth].max_bl_v);
    initBetaTable(beta_v, C.circuit[depth].max_bl_v, r_v[depth].begin(), F_ONE);
}

static ThreadSafeQueue<int> workerq,endq;
void pred1_worker(const vector<uniGate> &beg, vector<F>& uni_value,vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        for (auto gate=beg.begin()+l;gate!=beg.begin()+r;gate++) 
        {
            bool idx = gate->lu;
            uni_value[idx]+= beta_g[gate->g] * beta_u[gate->u] * gate->sc;
        }
        endq.Push(idx);
    }
}
void pred2_worker(const vector<binGate> &beg, vector<F>& bin_value,vector<F>& beta_g, vector<F>& beta_u,int*&L,int*&R) //F*& uni_value, layer &cur_layer,
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        for (auto gate=beg.begin()+l;gate!=beg.begin()+r;gate++) 
        {
            bin_value[gate->l] += beta_g[gate->g] * beta_u[gate->u] * beta_v[gate->v] * gate->sc;
        }
        endq.Push(idx);
    }
}

void verifier::predicatePhase1(u8 layer_id) 
{    
    auto &cur_layer = C.circuit[layer_id];
    uni_value[0].clear();
    uni_value[1].clear();
    vector<F>  univ[40];
    for(int i=0;i<32;i++)
    {
            univ[i].emplace_back(F(0));
            univ[i].emplace_back(F(0));
    }
    if(C.circuit[layer_id].uni_interval.size()>=2)
    {
        const int thd=32;
        int *L=new int [C.circuit[layer_id].uni_interval.size()],*R=new int [C.circuit[layer_id].uni_interval.size()];
        for (u64 j = 0; j <C.circuit[layer_id].uni_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur_layer.uni_interval[j].first;
                R[j]=cur_layer.uni_interval[j].second;
        }
            for(int i=0;i<thd;i++)
            { 
                thread t(pred1_worker, std::cref(cur_layer.uni_gates),std::ref(univ[i]),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R)); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(1));
            while(endq.Size()!=C.circuit[layer_id].uni_interval.size())
                this_thread::sleep_for (std::chrono::microseconds(1));
            endq.Clear();
            for(int i=0;i<thd;i++)
            {
                uni_value[0]+=univ[i][0];
                uni_value[1]+=univ[i][1];
            }
    }
    else for (auto &gate: cur_layer.uni_gates) 
    {
        bool idx = gate.lu;
        uni_value[idx]+= beta_g[gate.g] * beta_u[gate.u] * gate.sc;
    }
    bin_value[0] = bin_value[1] = bin_value[2] = F_ZERO;
}

void verifier::predicatePhase2(u8 layer_id) 
{
    uni_value[0] = uni_value[0] * beta_v[0];
    uni_value[1] = uni_value[1] * beta_v[0];

    auto &cur_layer = C.circuit[layer_id];
    vector<F>  binv[40];
    for(int i=0;i<32;i++)
    {
            binv[i].emplace_back(F(0));
            binv[i].emplace_back(F(0));
            binv[i].emplace_back(F(0));
    }
    if(C.circuit[layer_id].bin_interval.size()>=2)
    {
        const int thd=32;
        int *L=new int [C.circuit[layer_id].bin_interval.size()],*R=new int [C.circuit[layer_id].bin_interval.size()];
        for (u64 j = 0; j <C.circuit[layer_id].bin_interval.size(); ++j) 
        {
                workerq.Push(j);
                L[j]=cur_layer.bin_interval[j].first;
                R[j]=cur_layer.bin_interval[j].second;
        }
            for(int i=0;i<thd;i++)
            { 
                thread t(pred2_worker, std::cref(cur_layer.bin_gates),std::ref(binv[i]),std::ref(beta_g),std::ref(beta_u),std::ref(L),std::ref(R)); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(1));
            while(endq.Size()!=C.circuit[layer_id].bin_interval.size())
                this_thread::sleep_for (std::chrono::microseconds(1));
            endq.Clear();
            for(int i=0;i<thd;i++)
            {
                bin_value[0]+=binv[i][0];
                bin_value[1]+=binv[i][1];
                bin_value[2]+=binv[i][2];
            }
    }
    else for (auto &gate: cur_layer.bin_gates)
        bin_value[gate.l] += beta_g[gate.g] * beta_u[gate.u] * beta_v[gate.v] * gate.sc;
}

void verifier::prove(int commit_thread) 
{
    cout<<"Circuit Initiation finished"<<endl;
    cout<<"Proving Service started"<<endl;
    verifyGKR(); 
    verifyLasso(); 
    openCommit();
    cout<<"Our results:"<<endl;
    cout<<"Matrix multiplication Prover time: "<<matrix_time<<"s"<<endl;
    cout<<"Total Prover time: "<<prover_time + range_prover_time<<"s"<<endl;
    cout<<"Verifier time: "<<verifier_time<<"s"<<endl;
    cout<<"Proof size: "<<(double)(p->proof_size)/1024.0<<"KB"<<endl;
}

Fr * init_book_keeping(int m,int n,vector<Fr> &vec,int offset, vector<Fr> & ra)
{
    assert(ra.size()==m+n);
    Fr* pa[30];
    for(int j=0;j<=n;j++)
    {
        pa[j]=new Fr[1<<(m+n-j)];
    }
    int t=min((1<<(m+n)),(int)vec.size()-offset);
    memset(pa[0]+t,0,sizeof(Fr)*((1<<(m+n))-t));
    
    for(int j=0;j<t;j++)
        pa[0][j]=vec[j+offset];
    for(int j=1;j<=n;j++)
    {
        for(int kp=0;kp<(1<<(m+n-j));kp++)
        {
            pa[j][kp]=pa[j-1][kp+(1<<(m+n-j))]*ra[m+n-j]+pa[j-1][kp]*(1-ra[m+n-j]);
        }
    }
    for(int j=0;j<n;j++)
    {
        delete [] pa[j];
    }
    return pa[n];
}
Fr fm[1<<8][65536];
Fr * init_book_keeping_fast(int m,int n,int* vec,vector<Fr> & ra)
{
    assert(ra.size()==m+n);
    timer tt;
    
    Fr* ret=new Fr[1<<m];
    Fr *eqq=new Fr[1<<n],*eq1=new Fr[1<<(n/2)],*eq2=new Fr[1<<(n/2)];
    
    for(int j=0;j<(1<<n);j++)
    {
        eqq[j]=1;
        for(int i=0;i<n;i++)
        {
            if(j&(1<<(n-1-i)))
                eqq[j]*=ra[m+n-1-i];
            else
                eqq[j]*=1-ra[m+n-1-i];
        }
    }
    for(int j=0;j<(1<<(n/2));j++)
    {
        eq1[j]=1;
        for(int i=0;i<n/2;i++)
        {
            if(j&(1<<(n/2-1-i)))
                eq1[j]*=ra[m+n-1-i];
            else
                eq1[j]*=1-ra[m+n-1-i];
        }
    }
    for(int j=0;j<(1<<(n/2));j++)
    {
        eq2[j]=1;
        for(int i=0;i<n/2;i++)
        {
            if(j&(1<<(n/2-1-i)))
                eq2[j]*=ra[m+n/2-1-i];
            else
                eq2[j]*=1-ra[m+n/2-1-i];
        }
    }
    std::vector<std::thread> threads1,threads2;
    int num_threads=32;
    int chunk_size = (1<<(n/2))/num_threads;
    auto worker1 = [&](int start_id, int end_id) 
    {
        for(int R=start_id;R<end_id;R++)
        {
            fm[R][0]=0;
            fm[R][1]=eq2[R];
            for(int i=2;i<(1<<16);i++)
                fm[R][i]=fm[R][i-1]+eq2[R];
        }
    };
    for (int t = 0; t < num_threads; t++) 
    {
        int start_col = t * chunk_size;
        int end_col = min(start_col + chunk_size, (1<<(n/2)));
        threads1.emplace_back(worker1, start_col, end_col);
    }
    for (auto &th : threads1) 
    {
        if (th.joinable())
            th.join();
    }
    auto worker2 = [&](int start_id, int end_id) 
    {
        for(int kp=start_id;kp<end_id;kp++)
        {
            ret[kp]=0;
            Fr sum=0;
            for(int L=0;L<(1<<(n/2));L++)
            {
                Fr sumL=0;
                for(int R=0;R<(1<<(n/2));R++)
                {
                    int j=((L<<(n/2))+R)<<m;
                    if(vec[kp+j]!=0)
                        sumL+=fm[R][vec[kp+j]];
                }
                sum+=eq1[L]*sumL;
            }
            ret[kp]=sum;
        }
    };
    chunk_size = (1<<m)/num_threads;
    for (int t = 0; t < num_threads; t++) 
    {
        int start_col = t * chunk_size;
        int end_col = min(start_col + chunk_size, (1<<m));
        threads2.emplace_back(worker2, start_col, end_col);
    }
    for (auto &th : threads2) 
    {
        if (th.joinable())
            th.join();
    }
    return ret;
}

void sc_last_worker(Fr& A, Fr&B,Fr&C,Fr*& read_1,Fr*& read_2,Fr*& write_1,Fr*& write_2, int*& L,int*& R,Fr ran)
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        for(int kp=l;kp<r;kp++)
        {
            Fr & c=read_1[kp<<1|1], &d=read_1[kp<<1];
                if(c==d)
                    write_1[kp]=d;
                else
                    write_1[kp]=(c-d)*ran+d;
                Fr & a=read_2[kp<<1|1], &b=read_2[kp<<1];
                if(a==b)
                    write_2[kp]=b;
                else
                    write_2[kp]=(a-b)*ran+b;
                if(kp&1)
                {
                    bool az=write_2[kp-1].isZero(),bz=write_2[kp].isZero();
                    if(az && bz)
                        continue;
                    //if(!az)
                        A+=write_1[kp-1]*write_2[kp-1];
                    //if(!bz)
                        B+=write_1[kp]*write_2[kp];
                    C+=(write_1[kp]+write_1[kp]-write_1[kp-1])*(write_2[kp]+write_2[kp]-write_2[kp-1]);
                }
        }
        endq.Push(idx);
    }
}

void sc_last_worker_first(Fr& A, Fr&B,Fr&C,Fr*& read_1,Fr*& read_2, int*& L,int*& R)
{
    int idx;
    while (true)
    {
        bool ret=workerq.TryPop(idx);
        if(ret==false)
            return;
        int l=L[idx],r=R[idx];
        for(int k=l;k<r;k++)
        {
            if(k&1)
            {
                bool a_z=read_2[k].isZero(),b_z=read_2[k-1].isZero();
                if(a_z && b_z)
                    continue;
                B+=read_1[k]*read_2[k];
                C+=(read_1[k]+read_1[k]-read_1[k-1])*(read_2[k]+read_2[k]-read_2[k-1]); 
            }
            else if(!read_2[k].isZero())
                A+=read_1[k]*read_2[k];  
        }
        endq.Push(idx);
    }
}



pair<Fr,Fr> sum_check_product(Fr* f,Fr* g,int m,Fr* r,Fr ans)
{
    Fr A=0,B=0,C=0;
    Fr *zf=new Fr[1<<m], *zg=new Fr[1<<m];
    for(int k=0;k<(1<<m);k++)
    {
        zf[k]=2*f[k*2+1]-f[k*2];
        zg[k]=2*g[k*2+1]-g[k*2];
        C+=zf[k]*zg[k];
    }
    for(int k=0;k<(1<<(m+1));k++)
    {
        if(k&1)
            B+=f[k]*g[k];
        else
        {
            A+=f[k]*g[k];
        }
    }
    Fr a=(C+A)/2-B,b=2*B-C/2-3*A/2,c=A;
    Fr g1=a+b+c,g0=c;
    if(!(ans-(g1+g0)).isZero())
    {
        cout<<"ans!=g(0)+g(1) "<<m<<endl;
        exit(0);
    }
    if(!(a+b+c-B).isZero())
    {
        cout<<"sum check fail! order 1 "<<m<<endl;
        exit(0);
    }
    if(!(4*a+2*b+c-C).isZero())
    {
        cout<<"sum check fail! order 2 "<<m<<endl;
        exit(0);
    }
    Fr send_r=r[0];
    Fr * new_f=new Fr[1<<m], *new_g=new Fr[1<<m];
    for(int k=0;k<(1<<m);k++)
    {
        new_f[k]=(1-send_r)*f[k*2]+(send_r)*f[k*2+1];
        new_g[k]=(1-send_r)*g[k*2]+(send_r)*g[k*2+1];
    }
    Fr new_ans=send_r*send_r*a+send_r*b+c; //g(r)
    if(m==0)
    {
        assert((new_ans-new_f[0]*new_g[0]).isZero());
        return make_pair(new_f[0],new_g[0]);
    }
    else
        return sum_check_product(new_f,new_g,m-1,r+1,new_ans);
}


bool verifier::verifyGKR() 
{
    
    F alpha = F_ONE, beta = F_ZERO, relu_rou, final_claim_u1, final_claim_v1;
    r_u[C.size].resize(C.circuit[C.size - 1].bit_length);
    for (i8 i = 0; i < C.circuit[C.size - 1].bit_length; ++i)
        r_u[C.size][i].setByCSPRNG();
    vector<F>::const_iterator r_0 = r_u[C.size].begin();
    vector<F>::const_iterator r_1;

    auto previousSum = p->Vres(r_0, C.circuit[C.size - 1].size, C.circuit[C.size - 1].bit_length,C.size-1);
    timer ptimer,vtimer;
    timer mat_timer;
    ptimer.start();
    p -> sumcheckInitAll(r_0);
    ptimer.stop();
    prover_time+=ptimer.elapse_sec();
    int fc_cnt=0,now_fc_id;
    Fr prev_u1=0,prev_v1=0;
    for (int i = C.size - 1; i; --i)
    {
        if(C.circuit[i].ty==layerType::FCONN)
            fc_cnt++;
    } 
    
    timer init,sc;
    double init_t=0,sc_t=0;
    double matrix_tt_time=0;
    for (int i = C.size - 1; i; --i) 
    {
        timer tmp;
        tmp.start();
        auto &cur = C.circuit[i];
        ptimer.start();
        p->sumcheckInit(alpha, beta);
        ptimer.stop();
        prover_time+=ptimer.elapse_sec();
        
        // phase 1
        if (cur.zero_start_id < cur.size)
        {
            relu_rou.setByCSPRNG();
        }
        else 
            relu_rou = F_ONE;
        F previousRandom = F_ZERO;
        if(C.circuit[i].ty==layerType::FCONN)   //substitute matrix with thaler13
        {
            now_fc_id=--fc_cnt;
            int n,m,k;
            n=ceilPow2BitLength(p->fc_input_row[now_fc_id]);
            k=ceilPow2BitLength(p->fc_col[now_fc_id]);
            m=ceilPow2BitLength(p->fc_row[now_fc_id]);
            r_u[i].resize(m+n);
            r_v[i].resize(m+k);
            for (int j = m; j < m+n; ++j)
                r_u[i][j]=r_u[i+1][j-m+k];
            for (int j = m; j < m+k; ++j)
                r_v[i][j]=r_u[i+1][j-m];
            for(int j=0;j<m;j++)
            {
                r_u[i][j].setByCSPRNG();
                r_v[i][j]=r_u[i][j];
            }
        
            prev_u1=final_claim_u1;
            previousRandom = F_ZERO;
            Fr sum=0;
            ptimer.start();
            mat_timer.start();
            Fr* pa=init_book_keeping(m,n,p->val[0],p->fc_input_id[now_fc_id],r_u[i]);
            Fr* pb=init_book_keeping_fast(m,k,p->mat_val[now_fc_id],r_v[i]);
            pair<Fr,Fr> oracle=sum_check_product(pa,pb,m-1,r_u[i].data(),prev_u1);
            mat_timer.stop();
            ptimer.stop();
            matrix_time+=mat_timer.elapse_sec();
            prover_time+=ptimer.elapse_sec();    
            final_claim_u0[i]=oracle.first;
            final_claim_u1=0;
            final_claim_v0[i]=oracle.second;
            final_claim_v1=0;
        }
        else
        {
            timer normal_timer;
            normal_timer.start();
            ptimer.start();
            p->sumcheckInitPhase1(relu_rou);
            ptimer.stop();
            prover_time+=ptimer.elapse_sec();            
            r_u[i].resize(cur.max_bl_u);
            for (int j = 0; j < cur.max_bl_u; ++j) 
                r_u[i][j].setByCSPRNG();
            sc.start();
            for (i8 j = 0; j < cur.max_bl_u; ++j) 
            {
                F cur_claim, nxt_claim;
                ptimer.start();
                quadratic_poly poly = p->sumcheckUpdate1(previousRandom);
                ptimer.stop();
                prover_time+=ptimer.elapse_sec();     
                

                vtimer.start();
                cur_claim = poly.eval(F_ZERO) + poly.eval(F_ONE);
                nxt_claim = poly.eval(r_u[i][j]);

                if (cur_claim != previousSum) 
                {
                    cerr << cur_claim << ' ' << previousSum << endl;
                    fprintf(stderr, "Verification fail, phase1, circuit %d, current bit %d\n", i, j);
                    return false;
                }
                vtimer.stop();
                verifier_time+=vtimer.elapse_sec();
                previousRandom = r_u[i][j];
                previousSum = nxt_claim;
            }
            sc.stop();
            sc_t+=sc.elapse_sec();
            prev_u1=final_claim_u1;
            ptimer.start();
            p->sumcheckFinalize1(previousRandom, final_claim_u0[i], final_claim_u1);
            betaInitPhase1(i, alpha, beta, r_0, r_1, relu_rou);
            predicatePhase1(i);
            ptimer.stop();
            prover_time+=ptimer.elapse_sec();     
            
            if (cur.need_phase2) 
            {
                timer normal_timer2;
                normal_timer2.start();
                r_v[i].resize(cur.max_bl_v);
                for (int j = 0; j < cur.max_bl_v; ++j) 
                    r_v[i][j].setByCSPRNG();
              
                ptimer.start();
                p->sumcheckInitPhase2();
                ptimer.stop();
                prover_time+=ptimer.elapse_sec();     
                previousRandom = F_ZERO;
                sc.start();
                for (u32 j = 0; j < cur.max_bl_v; ++j) 
                {
                        ptimer.start();
                        quadratic_poly poly = p->sumcheckUpdate2(previousRandom);
                        ptimer.stop();
                        prover_time+=ptimer.elapse_sec();     
                        vtimer.start();
                        if (poly.eval(F_ZERO) + poly.eval(F_ONE) != previousSum) 
                        {
                            fprintf(stderr, "Verification fail, phase2, circuit level %d, current bit %d, total is %d\n", i, j,
                                    cur.max_bl_v);
                            return false;
                        }
                        vtimer.stop();
                        verifier_time+=vtimer.elapse_sec();
                        previousRandom = r_v[i][j];
                        previousSum = poly.eval(previousRandom);
                }
                sc.stop();
                sc_t+=sc.elapse_sec();
                prev_v1=final_claim_v1;
                ptimer.start();
                p->sumcheckFinalize2(previousRandom, final_claim_v0[i], final_claim_v1);
                betaInitPhase2(i);
                predicatePhase2(i);
                ptimer.stop();
                prover_time+=ptimer.elapse_sec();     
            }
            
        }
        
        
        vtimer.start();
        if (C.circuit[i].ty!=layerType::FCONN )  // for thaler13 no need to check this
        {
            F test_value = getFinalValue(final_claim_u0[i], final_claim_u1, final_claim_v0[i], final_claim_v1);
            if(previousSum != test_value)
            {
                std::cerr << test_value << ' ' << previousSum << std::endl;
                return false;
            }
        }
        
        

        if (~cur.bit_length_u[1])
            alpha.setByCSPRNG();
        else 
            alpha.clear();
        if ((~cur.bit_length_v[1]) || cur.ty == layerType::FFT)
            beta.setByCSPRNG();
        else 
            beta.clear();
        previousSum = alpha * final_claim_u1 + beta * final_claim_v1;
        vtimer.stop();
        verifier_time+=vtimer.elapse_sec();

        r_0 = r_u[i].begin();
        r_1 = r_v[i].begin();
        beta_u.clear();
        beta_v.clear();
        tmp.stop();
    }

    
    return true;
}

bool verifier::verifyLasso() 
{
    timer ptimer,vtimer;
    auto &cur = C.circuit[0];

    vector<F> sig_u(C.size - 1);
    for (int i = 0; i < C.size - 1; ++i) 
        sig_u[i].setByCSPRNG();
    vector<F> sig_v(C.size - 1);
    for (int i = 0; i < C.size - 1; ++i) 
        sig_v[i].setByCSPRNG();
    r_u[0].resize(cur.bit_length);
    for (int i = 0; i < cur.bit_length; ++i) 
        r_u[0][i].setByCSPRNG();
    auto r_0 = r_u[0].begin();
    
    F previousSum = F_ZERO;
    for (int i = 1; i < C.size; ++i) 
    {
        if (~C.circuit[i].bit_length_u[0])
            previousSum = previousSum + sig_u[i - 1] * final_claim_u0[i];
        if (~C.circuit[i].bit_length_v[0])
            previousSum = previousSum + sig_v[i - 1] * final_claim_v0[i];
    }
    ptimer.start();
    p->sumcheckLassoInit(sig_u, sig_v,r_u,r_v);
    ptimer.stop();
    prover_time+=ptimer.elapse_sec();
    F previousRandom = F_ZERO;
    int n=C.circuit[0].bit_length;
    Fr*pb1=new Fr[1<<(n+1)],* pa1[30],*pc1=pb1;
    Fr*pb2=new Fr[1<<(n+1)],* pa2[30],*pc2=pb2;
    for(int j=0;j<=n;j++)
    {
        pa1[j]=pc1;
        pc1+=1<<(n-j);
        pa2[j]=pc2;
        pc2+=1<<(n-j);
    }

    const __int128 FILTER=65535,NEG_FILTER=-65535;
    int *L=new int[1<<10],*R=new int[1<<10];
    ptimer.start();
    for(int i=0;i<=n;i++)
    {
        Fr A=0,B=0,C=0;
        if(i==0)
        {
            Fr* read_1=p->lasso_mult_v.data(),*read_2=p->val[0].data();
            const int k=10;
            const int thread_n=32; 
            int total_work=1<<n;
            int task_cnt=0;
            for (u64 j = 0; j < (1<<k); ++j) 
            {
                L[j]=(total_work>>k)*j;
                if(L[j]>=p->total_size[1])
                    break;
                R[j]=(total_work>>k)*(1+j);
                task_cnt++;
                workerq.Push(j);
            }
            Fr Aa[40],Bb[40],Cc[40];
            memset(Aa,0,sizeof(Aa));
            memset(Bb,0,sizeof(Bb));
            memset(Cc,0,sizeof(Cc));
            for(int j=0;j<thread_n;j++)
            {
                thread t(sc_last_worker_first,std::ref(Aa[j]),std::ref(Bb[j]),std::ref(Cc[j]),std::ref(read_1),std::ref(read_2),std::ref(L),std::ref(R)); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(1));
            while(endq.Size()!=task_cnt)
                this_thread::sleep_for (std::chrono::microseconds(1));
            endq.Clear();
            //Fr ax=0,bx=0,cx=0;
            for(int j=0;j<thread_n;j++)
            {
                A+=Aa[j];
                B+=Bb[j];
                C+=Cc[j];
            }
        }
        else if(i<=9)
        {
            timer tt;
            tt.start();
            Fr* read_1=p->lasso_mult_v.data(),*read_2=p->val[0].data(),*write_1=pa1[i],*write_2=pa2[i];
            if(i!=1)    
                read_1=pa1[i-1],read_2=pa2[i-1];
            const int k=10;
            const int thread_n=32; 
            int total_work=1<<(n-i);
            int taskcnt=0;
            for (u64 j = 0; j < (1<<k); ++j) 
            {
                L[j]=(total_work>>k)*j;
                R[j]=(total_work>>k)*(1+j);
                if((1<<i)*L[j]>=p->total_size[1])
                    break;
                ++taskcnt;
                workerq.Push(j);
            }
            Fr Aa[40],Bb[40],Cc[40];
            memset(Aa,0,sizeof(Aa));
            memset(Bb,0,sizeof(Bb));
            memset(Cc,0,sizeof(Cc));
            for(int j=0;j<thread_n;j++)
            {
                thread t(sc_last_worker,std::ref(Aa[j]),std::ref(Bb[j]),std::ref(Cc[j]),std::ref(read_1),std::ref(read_2),std::ref(write_1),std::ref(write_2),std::ref(L),std::ref(R),r_u[0][i-1]); 
                t.detach();
            }
            while(!workerq.Empty())
                this_thread::sleep_for (std::chrono::microseconds(10));
            while(endq.Size()!=taskcnt)
                this_thread::sleep_for (std::chrono::microseconds(10));
            endq.Clear();
            tt.stop();
            for(int j=0;j<thread_n;j++)
            {
                A+=Aa[j];
                B+=Bb[j];
                C+=Cc[j];
            }
        }
        else
        {
            Fr* read_1=pa1[i-1],*read_2=pa2[i-1],*write_1=pa1[i],*write_2=pa2[i];
            for(int kp=0;kp<(1<<(n-i));kp++)
            {
                Fr & c=read_1[kp<<1|1], &d=read_1[kp<<1];
                if(c==d)
                    write_1[kp]=d;
                else
                    write_1[kp]=(c-d)*r_u[0][i-1]+d;
                Fr & a=read_2[kp<<1|1], &b=read_2[kp<<1];
                if(a==b)
                    write_2[kp]=b;
                else
                    write_2[kp]=(a-b)*r_u[0][i-1]+b;
                if(kp&1)
                {
                    B+=write_1[kp]*write_2[kp];
                    C+=(write_1[kp]+write_1[kp]-write_1[kp-1])*(write_2[kp]+write_2[kp]-write_2[kp-1]);
                }
                else
                    A+=write_1[kp]*write_2[kp];
            }
        }
        vtimer.start();
        Fr aa=(C+A)/2-B,bb=2*B-C/2-3*A/2,cc=A;
        Fr g1=aa+bb+cc,g0=cc;
        if(i!=n)
            assert(previousSum==g0+g1);
        if(i!=n)
            previousSum=aa*r_u[0][i]*r_u[0][i]+bb*r_u[0][i]+cc; 
        vtimer.stop();
        verifier_time+=vtimer.elapse_sec();
    }
    
    ptimer.stop();
    prover_time+=ptimer.elapse_sec();
    F gr = F_ZERO;

    eval_in=pa2[n][0];
    
    beta_g.resize(1ULL << cur.bit_length);
    ptimer.start();
    
    initBetaTable(beta_g, cur.bit_length, r_0, F_ONE);
    for (int i = 1; i < C.size; ++i) 
    {
        if (~C.circuit[i].bit_length_u[0]) 
        {
            beta_u.resize(1ULL << C.circuit[i].bit_length_u[0]);
            initBetaTable(beta_u, C.circuit[i].bit_length_u[0], r_u[i].begin(), sig_u[i - 1]);
            for (u32 j = 0; j < C.circuit[i].size_u[0]; ++j)
                gr = gr + beta_g[C.circuit[i].ori_id_u[j]] * beta_u[j];
        }

        if (~C.circuit[i].bit_length_v[0]) 
        {
            beta_v.resize(1ULL << C.circuit[i].bit_length_v[0]);
            initBetaTable(beta_v, C.circuit[i].bit_length_v[0], r_v[i].begin(), sig_v[i - 1]);
            for (u32 j = 0; j < C.circuit[i].size_v[0]; ++j)
                gr = gr + beta_g[C.circuit[i].ori_id_v[j]] * beta_v[j];
        }
    }
    ptimer.stop();
    
    prover_time+=ptimer.elapse_sec();
    
    beta_u.clear();
    beta_v.clear();

    vtimer.start();
    assert (eval_in * gr == previousSum);
    vtimer.stop();
    verifier_time+=vtimer.elapse_sec();

    fprintf(stderr, "All verification passed!!\n");

    beta_g.clear();
    beta_gs.clear();
    beta_u.clear();
    beta_v.clear();
    r_u.resize(1);
    r_v.clear();

    sig_u.clear();
    sig_v.clear();
    return true;
}
__int128 convertx(Fr x)	
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
    __int128 V=0;	
    for(int j=size-1;j>=0;j--)
    {	
        V=(V<<8)+bf[j];	
    }
    if(sign)	
        V=-V;	
    return V;	
}
bool verifier::openCommit() 
{
    Fr*Lp=new Fr[1<<(comm.l/2)],*Rp=new Fr[1<<(comm.l-comm.l/2)];
    timer verify_input;
    verify_input.start();
    brute_force_compute_LR(Lp,Rp,r_u[0].data(),comm.l);
    int rownum=(1<<comm.l/2),colnum=(1<<(comm.l-comm.l/2));
    
    ll* row=new ll[1<<comm.l];
    pair<double,double> tim=hyrax::open(comm.w,r_u[0].data(),eval_in,comm.G,comm.g,Lp,Rp,comm.comm,comm.l);
    prover_time+=tim.first;
    verifier_time+=tim.second;
    return true;
}
