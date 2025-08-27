//
// Created by 69029 on 3/16/2021.
//
#undef NDEBUG
#include "neuralNetwork.hpp"
#include "utils.hpp"
#include "global_var.hpp"
#include <polynomial.h>
#include <circuit.h>
#include <iostream>
#include <cmath>

using std::cerr;
using std::endl;
using std::max;
using std::ifstream;
using std::ofstream;

ifstream in;
ifstream conf;
ofstream out;
namespace multi_max 
{
template<class T>
T max(T head) {
    return head;
}
template<class T, typename... Args>
T max(T head, Args... args) {
    T t = max<T>(args...);
    return (head > t)?head:t;
}
} // end of namespace

neuralNetwork::neuralNetwork(i64 psize_x, i64 psize_y, i64 pchannel, i64 pparallel, const string &i_filename,
                             const string &c_filename, const string &o_filename, bool i_llm,
                             int headnum, int headdim, int attn_dim, int linear_dim) :
        pic_size_x(psize_x), pic_size_y(psize_y), pic_channel(pchannel), pic_parallel(pparallel),
        SIZE(0), NCONV_FAST_SIZE(1), NCONV_SIZE(2), FFT_SIZE(5),
        AVE_POOL_SIZE(1), FC_SIZE(1), RELU_SIZE(2), act_ty(RELU_ACT),
        headnum(headnum), headdim(headdim), attn_dim(attn_dim), linear_dim(linear_dim) 
    {
        is_llm=i_llm;
        in.open(i_filename);
        conf.open(c_filename);
}


void Out_group(Fr a,Fr b,Fr c)
{ 
    char flaga=' ',flagb=' ',flagc=' ';
    if (a>1e10)
    {
        flaga='-';
        a=-a;
    }
    
    if (b>1e10)
    {
        b=-b;
        flagb='-';
    }
    if (c>1e10)
    {
        flagc='-';
        c=-c;
    }
}

// input:   [data]
//          [[conv_kernel || relu_conv_bit_decmp]{sec.size()}[max_pool]{if maxPool}[pool_bit_decmp]]{conv_section.size()}
//          [fc_kernel || relu_fc_bit_decmp]
void neuralNetwork::initParam(prover &pr,int depth) 
{
    total_in_size = 0;
    total_para_size = 0;
    total_relu_in_size = 0;
    total_ave_in_size = 0;
    total_max_in_size = 0;
    // data
    i64 pos = 32*1024;
    for (int i = 0; i < full_conn.size(); ++i) 
    {
        auto &fc = full_conn[i];
        refreshFCParam(fc);
        // fc_kernel
        pr.fc_real_input_row[i]=len;
        len=pr.fc_input_row[i]=1<<ceilPow2BitLength(len);
        pr.fc_real_input_col[i]=fc.channel_in;
        pr.fc_real_row[i]=fc.channel_in;
        pr.fc_input_col[i]=1<<ceilPow2BitLength(fc.channel_in);
        
        fc.channel_in=pr.fc_row[i]=1<<ceilPow2BitLength(fc.channel_in);
        pr.fc_real_col[i]=fc.channel_out;
        fc.channel_out=pr.fc_col[i]=1<<ceilPow2BitLength(fc.channel_out);
        
        fc.weight_start_id = pos;   // TODO calc FC pos
        pr.fc_start_id[i]=fc.weight_start_id;
        u32 para_size = pr.fc_row[i] * pr.fc_col[i];
        pos += para_size;
        total_para_size += para_size;
        //fc.bias_start_id = pos;
        //pos += channel_out;
        //total_para_size += channel_out;
        //fprintf(stderr, "full conn  bias   weight: %11lld%11lld\n", channel_out, total_para_size);
    }
    total_in_size = pos;
    vector<string> layers={"l1","l2","l3","fcon","round","MHA_QK","softmax*V","softmax*v","soft_end","fcon2","round2","l1","l2","l3","fcon3","round","gelu1","g2","g3","fcon4","round"};
    SIZE=1+layers.size()*depth ;         //TODO here is very important to avoid strange memory errors
}
void neuralNetwork::merge_layer(prover &pr,i64 layer_id)
{
    int cntp=0;
    for(int i=0;i<layer_id;i++)
    {
        if( (int)pr.C.circuit[i].ty==4 || (int)pr.C.circuit[i].ty==11)
        {
            ++cntp;
            layer v;
            v.ty=pr.C.circuit[i].ty;
            pr.C.circuit.push_back(v);
            swap(pr.C.circuit[i],pr.C.circuit[pr.C.circuit.size()-1]);
            vector<F> f;
            pr.val.push_back(f);
            swap(pr.val[i],pr.val[pr.val.size()-1]);
        }
    }
    vector<vector<F> >::iterator itf=pr.val.begin(); 
    int cnt2=0;
    for (vector<layer>::iterator it = pr.C.circuit.begin(); it <pr.C.circuit.end() ;) 
    {
        if ((int)it->ty==4 || (int)it->ty==11) 
        {
            if (cnt2<cntp)
            {
                it = pr.C.circuit.erase(it);
                itf=pr.val.erase(itf);
                cnt2++;
            }
            else
            {
                ++it;
                ++itf;
            }
        } 
        else 
        {
            ++it;
            ++itf;
        }
    }
    int offset[5]={0,0,0,0,0};
    int lazy_offset[5]={0,0,0,0,0};
    for(int i=4;i<layer_id;i++)
    {
        auto t=pr.C.circuit[i].ty;
        int os=0;
        if(t==layerType::MHA_QK || t==layerType::GELU_1 || t==layerType::LAYER_NORM_1)
            os=1;
        else if(t==layerType::SOFTMAX_1 || t==layerType::GELU_2 || t==layerType::LAYER_NORM_2)
            os=2;
        else if(t==layerType::SOFTMAX_2 || t==layerType::GELU_3 || t==layerType::LAYER_NORM_3)
            os=3;
        else if(t==layerType::SOFTMAX_3)
            os=4;
        else
            break;
        assert(os!=0);
        assert(pr.val[i].size()==pr.C.circuit[i].size);
        
        for(int j=0;j<pr.C.circuit[i].size;j++)
        {
            pr.val[os].emplace_back(pr.val[i][j]);
        }
    }
    for(int i=1;i<4;i++)
    {
        auto t=pr.C.circuit[i].ty;
        int os=0, os1=0;
        if(t==layerType::MHA_QK || t==layerType::GELU_1 || t==layerType::LAYER_NORM_1)
            os=1;
        else if(t==layerType::SOFTMAX_1 || t==layerType::GELU_2 || t==layerType::LAYER_NORM_2)
            os=2;
        else if(t==layerType::SOFTMAX_2 || t==layerType::GELU_3 || t==layerType::LAYER_NORM_3)
            os=3;
        else if(t==layerType::SOFTMAX_3)
            os=4;
        if(os==i)
        {
            pr.C.circuit[os].uni_interval.emplace_back(make_pair(0,pr.C.circuit[i].uni_gates.size()));
            pr.C.circuit[os].bin_interval.emplace_back(make_pair(0,pr.C.circuit[i].bin_gates.size()));
        }
    }
    for(int i=4;i<layer_id;i++)
    {
        auto t=pr.C.circuit[i].ty;
        int os=0, os1=0;
        if(t==layerType::MHA_QK || t==layerType::GELU_1 || t==layerType::LAYER_NORM_1)
            os=1;
        else if(t==layerType::SOFTMAX_1 || t==layerType::GELU_2 || t==layerType::LAYER_NORM_2)
            os=2;
        else if(t==layerType::SOFTMAX_2 || t==layerType::GELU_3 || t==layerType::LAYER_NORM_3)
            os=3;
        else if(t==layerType::SOFTMAX_3)
            os=4;
        else
            break;
        assert(os!=0);
        
        auto t2=pr.C.circuit[i-1].ty;
        if(t2==layerType::MHA_QK || t2==layerType::GELU_1 || t2==layerType::LAYER_NORM_1)
            os1=1;
        else if(t2==layerType::SOFTMAX_1 || t2==layerType::GELU_2 || t2==layerType::LAYER_NORM_2)
            os1=2;
        else if(t2==layerType::SOFTMAX_2 || t2==layerType::GELU_3 || t2==layerType::LAYER_NORM_3)
            os1=3;
        else if (t2==layerType::SOFTMAX_3)
            os1=4;
        
        assert(i!=os);
        offset[os1]=offset[os]=0;
        for(int j=1;j<i;j++)
        {
            auto tp=pr.C.circuit[j].ty;
            int tos=0;
            if(tp==layerType::MHA_QK || tp==layerType::GELU_1 || tp==layerType::LAYER_NORM_1)
                tos=1;
            else if(tp==layerType::SOFTMAX_1 || tp==layerType::GELU_2 || tp==layerType::LAYER_NORM_2)
                tos=2;
            else if(tp==layerType::SOFTMAX_2 || tp==layerType::GELU_3 || tp==layerType::LAYER_NORM_3)
                tos=3;
            else if(tp==layerType::SOFTMAX_3)
                tos=4;
            if(tos==os)
                offset[os]+=pr.C.circuit[j].size;
        }
        for(int j=1;j<i-1;j++)
        {
            auto tp=pr.C.circuit[j].ty;
            int tos=0;
            if(tp==layerType::MHA_QK || tp==layerType::GELU_1 || tp==layerType::LAYER_NORM_1)
                tos=1;
            else if(tp==layerType::SOFTMAX_1 || tp==layerType::GELU_2 || tp==layerType::LAYER_NORM_2)
                tos=2;
            else if(tp==layerType::SOFTMAX_2 || tp==layerType::GELU_3 || tp==layerType::LAYER_NORM_3)
                tos=3;
            else if(tp==layerType::SOFTMAX_3)
                tos=4;
            if(tos==os1)
                offset[os1]+=pr.C.circuit[j].size;
        }
        pr.C.circuit[os].uni_interval.emplace_back(make_pair(pr.C.circuit[os].uni_gates.size(),pr.C.circuit[os].uni_gates.size()+pr.C.circuit[i].uni_gates.size()));
        for(auto g=pr.C.circuit[i].uni_gates.begin();g!=pr.C.circuit[i].uni_gates.end();g++)
        {
            if((int)g->lu==0)
                pr.C.circuit[os].uni_gates.emplace_back(g->g+offset[os], g->u,0,g->sc);
            else
            {
                pr.C.circuit[os].uni_gates.emplace_back(g->g+offset[os], g->u+offset[os1],os1,g->sc);
                assert(os==os1+1);
                assert(g->u<pr.C.circuit[i-1].size);
                assert(g->u+offset[os1]<pr.val[os1].size());
            }
        }
        pr.C.circuit[os].bin_interval.emplace_back(make_pair(pr.C.circuit[os].bin_gates.size(),pr.C.circuit[os].bin_gates.size()+pr.C.circuit[i].bin_gates.size()));
        for(auto g=pr.C.circuit[i].bin_gates.begin();g!=pr.C.circuit[i].bin_gates.end();g++)
        {
            if((int)g->l==0)
                pr.C.circuit[os].bin_gates.emplace_back(g->g+offset[os],g->u,g->v,g->sc,g->l);
            else if((int)g->l==1)
                pr.C.circuit[os].bin_gates.emplace_back(g->g+offset[os],g->u+offset[os1],g->v+offset[os1],g->sc,g->l);
            else if((int)g->l==2)
                pr.C.circuit[os].bin_gates.emplace_back(g->g+offset[os],g->u+offset[os1],g->v,g->sc,g->l);
            assert(g->g+offset[os]<pr.val[os].size());
        }
    }
    pr.C.circuit.erase(pr.C.circuit.begin()+5,pr.C.circuit.begin()+layer_id-cntp);
    pr.val.erase(pr.val.begin()+5,pr.val.begin()+layer_id-cntp);

    
    pr.C.size=pr.C.circuit.size();
    for(int i=1;i<pr.C.size;i++)
    {
        initLayer(pr.C.circuit[i], pr.val[i].size(), pr.C.circuit[i].ty);
        if(pr.C.circuit[i].ty!=layerType::FCONN)
            checkNormalLayer(pr.C.circuit[i],i,pr.val);
    }
}
void neuralNetwork::create(prover &pr, bool merge) 
{
    compute_e_table();
    initParam(pr,layer_num);
    pr.C.init(Q_BIT_SIZE, SIZE);
    pr.val.resize(SIZE);
    val = pr.val.begin();
    i64 layer_id = 0;
    inputLayer(pr.C.circuit[layer_id++]);
    for (int i = 0; i < full_conn.size(); ++i) 
    {
        auto &fc = full_conn[i];
        refreshFCParam(fc);
        readFconWeight(fc.weight_start_id,pr.fc_real_row[i],pr.fc_real_col[i],i);
    }
    timer T;
    int logn = pr.C.circuit[0].bit_length;
    u64 n_sqrt = 1ULL << (logn - (logn >> 1));
    int c=0;
    pr.gens.resize(n_sqrt);
    G1 base=gen_gi(pr.gens.data(),n_sqrt);
    pr.gens.push_back(base);
    T.start();
    pr.commitInput(pr.gens,32);  //commit weight
    T.stop();
    pr.proof_size+= 1<<(pr.cc.l/2);
    cout<<"Model weight commit time: "<<T.elapse_sec()<<"s"<<endl;
    cout<<"Start initiating circuit"<<endl;
    for (int i = 0; i < full_conn.size(); ++i) 
    {
        auto &fc = full_conn[i];
        refreshFCParam(fc);
        if(i==0)
        {
            q_offset=0;
        }
        read_layer_norm(0);
        bool * sparsity=new bool[pr.fc_input_row[i]*pr.fc_input_col[i]];
        memset(sparsity,0,sizeof(bool)*pr.fc_input_row[i]*pr.fc_input_col[i]);
        int cnt=0;
        for(int j=0;j<pr.fc_real_input_row[i];j++)
        for(int k=0;k<pr.fc_real_input_col[i];k++)
        {
            ++cnt;
            sparsity[j*pr.fc_input_col[i]+k]=true;
        }
        if(i%4==0 || i%4==2)
        {
            ln_checker_layer1(pr.C.circuit[layer_id], layer_id,0,input_e,input_c,pr.fc_real_input_col[0],sparsity);
            ln_checker_layer2(pr.C.circuit[layer_id], layer_id,0,input_e,input_c,sparsity);
            ln_checker_layer3(pr.C.circuit[layer_id], layer_id,0,input_e,input_c,pr.fc_real_input_col[0],sparsity);
        }
        pr.fc_input_id[i]=q_offset;
        fullyConnLayer(pr.C.circuit[layer_id], layer_id, fc.weight_start_id,q_offset,0);
        int cx=7,ex=-8,cy=3,ey=-8;
        float c_A,e_A,c_B=1,e_B=-10,c_C=7,e_C=-8;
        c_A=input_c;
        e_A=input_e;
        
        roundLayer(pr.C.circuit[layer_id], layer_id,(float)c_A*c_B/c_C*pow(2,e_A+e_B-e_C));
        if(i%4==0)
        {
            multi_head_matrix_QK(pr.C.circuit[layer_id], layer_id);
            softmax_layer_1(pr.C.circuit[layer_id], layer_id,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(1,-8));
            softmax_layer_2(pr.C.circuit[layer_id], layer_id,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(1,-8));
            softmax_layer_3(pr.C.circuit[layer_id], layer_id,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(2,e_C)*c_C,pow(1,-8));
        }
        if(i%4==2)
        {
            gelu_checker_layer1(pr.C.circuit[layer_id], layer_id,pr.fc_real_col[i],-8,48, -8,217,-8, 252,-8, 615,ex,cx,ey,cy);
            gelu_checker_layer2(pr.C.circuit[layer_id], layer_id,pr.fc_real_col[i],-8,48, -8,217,-8, 252,-8, 615,ex,cx,ey,cy);
            gelu_checker_layer3(pr.C.circuit[layer_id], layer_id,pr.fc_real_col[i],-8,48, -8,217,-8, 252,-8, 615,ex,cx,ey,cy);
        }
    }
   
    if(merge)
    {
        merge_layer(pr,layer_id);
    }
    
    total_in_size += total_max_in_size + total_ave_in_size + total_relu_in_size;
    initLayer(pr.C.circuit[0], total_in_size, layerType::INPUT);
    assert(total_in_size == pr.val[0].size());

    pr.C.initSubset();
    
    int cnt=0;
    for(int i=0;i<pr.C.size;i++)
        cnt+=pr.val[i].size();
    int bin=0,uni=0;
    for(int i=0;i<pr.C.size;i++)
    {
        bin+=pr.C.circuit[i].bin_gates.size();
        uni+=pr.C.circuit[i].uni_gates.size();
    }
    pr.mat_val=mat_values;
}

void neuralNetwork::inputLayer(layer &circuit) 
{
    initLayer(circuit, total_in_size, layerType::INPUT);
    for (i64 i = 0; i < total_in_size; ++i) 
        circuit.uni_gates.emplace_back(i, 0, 0, 1);

    calcInputLayer(circuit);
}


pair<int,int> search(double scale)
{
    double mindiff=1e9;
    int best_e,best_c;
    for(int e=-10;e<=10;e++)
    for(int c=1;c<=800;c++)
    {
        double s=pow(2,e)*c;
        if(abs(s-scale)<mindiff)
        {
            mindiff=abs(s-scale);
            best_e=e;
            best_c=c;
        }
    }
    return make_pair(best_e,best_c);
}

void neuralNetwork::read_layer_norm(int ln_id)
{
    //int layer_norm_w_c[30], layer_norm_w_e[30],layer_norm_b_c[30], layer_norm_b_e[30];
    //int layer_norm_w_q_start[30],layer_norm_b_q_start[30];
    int orgsize=val[0].size() ;
    val[0].resize(orgsize+2*channel_out);
    total_in_size +=2*channel_out;
    layer_norm_w_c[ln_id]=1;
    layer_norm_w_e[ln_id]=0;
    layer_norm_b_c[ln_id]=1;
    layer_norm_b_e[ln_id]=-8;
    layer_norm_w_q_start[ln_id]=orgsize;
    layer_norm_b_q_start[ln_id]=orgsize+channel_out;

    for(int i=0;i<channel_out;i++)
        val[0][i+layer_norm_w_q_start[ln_id]]=1;
    for(int i=0;i<channel_out;i++)
        val[0][i+layer_norm_b_q_start[ln_id]]=1;

}   
std::ostream& operator<<(std::ostream& os, __int128_t value) {
    if (value < 0) {
        os << '-';
        value = -value;
    }
    // save flags to restore them
    std::ios_base::fmtflags flags(os.flags());
    // set zero fill
    os << std::setfill('0') << std::setw(13);

    // 128-bit number has at most 39 digits,
    // so the below loop will run at most 3 times
    const int64_t modulus = 10000000000000; // 10**13
    do {
        int64_t val = value % modulus;
        value /= modulus;
        if (value == 0) {
            os.flags(flags);
            return os << val;
        }
        os << val;
    } while (1);
}
void neuralNetwork::ln_checker_layer1(layer &circuit, i64 &layer_id, int ln_id, int ey,int cy,int real_cn_in,bool* sparsity_map)
{
    
    int cw= layer_norm_w_c[ln_id]; //scale w  //TODO we need to fix all layer_norm value's read
    int ew= layer_norm_w_e[ln_id];
    int cb= layer_norm_b_c[ln_id];  //scale b
    int eb= layer_norm_b_e[ln_id];
    double sw=pow(2,ew)*cw,sb=pow(2,eb)*cb;
    double sy=pow(2,ey)*cy;
    int qw_off= layer_norm_w_q_start[ln_id];  // place of w vector
    int qb_off= layer_norm_b_q_start[ln_id];  // place of b vector
    
    int c1,e1,c2,e2;
    pair<int,int> S1,S2;
    S1=search(sw*sqrt(real_cn_in)/sy); // TODO: s1 is wrong, channel_out should be something else
    S2=search(sb/sy);
    e1=S1.first;
    c1=S1.second;
    e2=S2.first;
    c2=S2.second;
    layer_norm_c1[ln_id]=c1;
    layer_norm_e1[ln_id]=e1;
    layer_norm_c2[ln_id]=c2;
    layer_norm_e2[ln_id]=e2;
    int m=multi_max::max(1,-e1,-e2);
    i64 block_len = len* channel_in;
    int output_size=block_len+3*len;
    initLayer(circuit, output_size, layerType::LAYER_NORM_1);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=block_len+2*len;
    int orgsize=val[0].size() ;
    val[0].resize(orgsize+10+block_len*2+4*len);// y(768*len), sum,B, sigma, delta1, delta2(768*len)
    for(int i=orgsize;i<val[0].size();i++)
        val[0][i]=0;
    ln_aux_start=orgsize; 
    val[0][orgsize]=1; 
    total_relu_in_size += 10+ block_len*2+4*len;
    ll qx[1024]={0},a[1024]={0};
    ll qw[1024],qb[1024];
    for (i64 co = 0; co < channel_in; ++co) 
    {
        qw[co]=convert(val[0][qw_off+co]);
        qb[co]=convert(val[0][qb_off+co]);
    }
    __int128 mn=1;
    for (i64 i = 0; i < len; i++)
    {
        ll sum=0,B=1,sigma;
        if(!sparsity_map[i*channel_in])
        {
            continue;
        }
        for (i64 co = 0; co < channel_in; ++co) 
        {
            i64 g = matIdx(i, co, channel_in);
            qx[co]=convert(val[0][q_offset+g]);
            sum+=qx[co];
        }
        for (i64 co = 0; co < channel_in; ++co) 
        {
            if(!sparsity_map[i*channel_in+co])
                continue;
            a[co]=qx[co]*real_cn_in-sum;  //TODO here has to change
            B+=a[co]*a[co];
        }
        assert(B!=0);
        sigma=round(sqrt(B));
        assert(sigma!=0);
        assert(sigma*sigma+sigma+1-B<(1ll<<32) && sigma*sigma+sigma+1-B>0);
        assert(B-sigma*sigma+sigma<(1ll<<32) && B-sigma*sigma+sigma>0);
        Fr delta1= Fr(sigma*sigma+sigma+1-B)*Fr(B-sigma*sigma+sigma);
        for (i64 co = 0; co < channel_in; ++co) 
        {
            i64 g = matIdx(i, co, channel_in);
            if(!sparsity_map[g])
            {
                continue;
            }
            ll qy = round(pow(2,e1)*c1*qw[co]*a[co]/sigma+pow(2,e2)*c2*qb[co]);
            
            int y_off=orgsize+10+g;
            int d2_off=orgsize+10+block_len+len*4+matIdx(i, co, real_cn_in);
            val[0][y_off]=qy;
            ll term1,term2;
            term1=(2*qy+1)*(1ll<<(m-1))*sigma+1-(1ll<<(e1+m))*c1*qw[co]*a[co]-(1ll<<(e2+m))*c2*qb[co]*sigma;
            term2=(1ll<<(e1+m))*c1*qw[co]*a[co]+(1ll<<(e2+m))*c2*qb[co]*sigma-(2*qy-1)*(1ll<<(m-1))*sigma+1;
            assert(term1>0&&term2>0);
            val[0][d2_off]=Fr(term1)*Fr(term2);
            positive_check+=1;  //add one d2
        }
        int sum_off=orgsize+10+block_len+i;
        val[0][sum_off]=sum;
        int b_off=orgsize+10+block_len+len+i;
        val[0][b_off]=B;
        int sig_off=orgsize+10+block_len+len*2+i;
        val[0][sig_off]=sigma;
        int d1_off=orgsize+10+block_len+len*3+i;
        val[0][d1_off]=delta1;
        positive_check+=1;  //add one d1
    }
    Fr SUM=0;
    for (i64 i = 0; i < len; i++)
    {
        if(!sparsity_map[i*channel_in])
        {
            continue;
        }
        int sum_off=orgsize+10+block_len+i;
        int sig_off=orgsize+10+block_len+len*2+i;
        int b_off=orgsize+10+block_len+len+i;
        for (i64 co = 0; co < channel_in; ++co) 
        {
            int g = matIdx(i, co, channel_in);
            if(!sparsity_map[g])
            {
                continue;
            }
            circuit.uni_gates.emplace_back(g, q_offset+g, 0, real_cn_in); // verify a with x and sum 
            circuit.uni_gates.emplace_back(g, sum_off, 0, -1);
            if(i==2)
            {
                SUM+=val[0][q_offset+g];
            }
            circuit.uni_gates.emplace_back(block_len+len*2+i, q_offset+g, 0, 1); //+xi
        }
        circuit.bin_gates.emplace_back(block_len+i, sig_off,sig_off, 1,0);//sigma^2
        circuit.uni_gates.emplace_back(block_len+i, sig_off, 0, 1);  //sigma
        circuit.uni_gates.emplace_back(block_len+i, orgsize, 0, 1); //+1
        circuit.uni_gates.emplace_back(block_len+i, b_off, 0, -1); //-B


        circuit.bin_gates.emplace_back(block_len+len+i, sig_off,sig_off, -1,0);//-sigma^2
        circuit.uni_gates.emplace_back(block_len+len+i, sig_off, 0, 1);  //sigma
        circuit.uni_gates.emplace_back(block_len+len+i, b_off, 0, 1); //+B
        circuit.uni_gates.emplace_back(block_len+len*2+i, sum_off, 0, -1);  //-SUM
    }
    
    calcNormalLayer(circuit, layer_id);
    for (i64 i = 0; i < len; i++)
    {
        assert(val[layer_id][block_len+len*2+i].isZero());
    }

    layer_id++;
}

void neuralNetwork::ln_checker_layer2(layer &circuit, i64 &layer_id, int ln_id, int ey,int cy,bool* sparsity_map)
{
    int cw= layer_norm_w_c[ln_id]; //scale w  //TODO we need to fix all layer_norm value's read
    int ew= layer_norm_w_e[ln_id];
    int cb= layer_norm_b_c[ln_id];  //scale b
    int eb= layer_norm_b_e[ln_id];
    int qw_off= layer_norm_w_q_start[ln_id];  // place of w vector
    int qb_off= layer_norm_b_q_start[ln_id];  // place of b vector
    int c1,e1,c2,e2;
    c1=layer_norm_c1[ln_id];
    e1=layer_norm_e1[ln_id];
    c2=layer_norm_c2[ln_id];
    e2=layer_norm_e2[ln_id];
    int m=multi_max::max(-e1,-e2,1);
    i64 block_len = len* channel_in;
    int output_size=2*block_len+2*len;
    initLayer(circuit, output_size, layerType::LAYER_NORM_2);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=2*block_len;
    int orgsize=ln_aux_start;
    //ll qx[1024],a[1024];
    //ll qw[1024],qb[1024];
    for (i64 i = 0; i < len; i++)
    {
        int sum_off=orgsize+10+block_len+i;
        int sig_off=orgsize+10+block_len+len*2+i;
        int b_off=orgsize+10+block_len+len+i;

        // (2*qy+1)*(1ll<<(m-1))*sigma+1  -(1ll<<(e1+m))*c1*qw[co]*a[co]  -(1ll<<(e2+m))*c2*qb[co]*sigma;

        ll SUM=0,cir=0;
        for (i64 co = 0; co < channel_in; ++co) 
        {
            int g = matIdx(i, co, channel_in);
            int y_off=orgsize+10+g;
            if(!sparsity_map[g])
            {
                continue;
            }
            circuit.bin_gates.emplace_back(g, g, qw_off+co, -(1ll<<(e1+m))*c1 ,2 ); // -qw[co]*a[co]
            circuit.bin_gates.emplace_back(g, sig_off, qb_off+co, -(1ll<<(e2+m))*c2 ,0); //-(1ll<<(e2+m))*c2*qb[co]*sigma;
            circuit.bin_gates.emplace_back(g, sig_off,y_off ,1<<m ,0 );
            circuit.uni_gates.emplace_back(g, orgsize,0 ,1 );
            circuit.uni_gates.emplace_back(g, sig_off,0 ,(1<<(m-1)) );

            circuit.bin_gates.emplace_back(g+block_len, g, qw_off+co, (1ll<<(e1+m))*c1 ,2 ); // qw[co]*a[co]
            circuit.bin_gates.emplace_back(g+block_len, sig_off, qb_off+co, (1ll<<(e2+m))*c2 ,0); //(1ll<<(e2+m))*c2*qb[co]*sigma;
            circuit.bin_gates.emplace_back(g+block_len, sig_off,y_off ,-(1<<m) ,0 );
            circuit.uni_gates.emplace_back(g+block_len, orgsize,0 ,1 );
            circuit.uni_gates.emplace_back(g+block_len, sig_off,0 ,(1<<(m-1)) );

            circuit.bin_gates.emplace_back(2*block_len+len+i, g, g ,1,1 );
            cir=convert(val[layer_id-1][g]);
            SUM+=cir*cir;
        }
        if(!sparsity_map[i*channel_in])
        {
            continue;
        }
        ll P=convert(val[0][b_off]);
        circuit.bin_gates.emplace_back(2*block_len+i, block_len+i,block_len+len+i, 1,1);
        int d1_off=orgsize+10+block_len+len*3+i;
        circuit.uni_gates.emplace_back(2*block_len+i, d1_off, 0, -1); 
        circuit.uni_gates.emplace_back(2*block_len+len+i, orgsize, 0, 1); 
        circuit.uni_gates.emplace_back(2*block_len+len+i, b_off, 0, -1);  //b=sigma^2+1
        /*circuit.bin_gates.emplace_back(block_len+i, sig_off,sig_off, 1,0);//sigma^2
        circuit.uni_gates.emplace_back(block_len+i, sig_off, 0, 1);  //sigma
        circuit.uni_gates.emplace_back(block_len+i, orgsize, 0, 1); //+1
        circuit.uni_gates.emplace_back(block_len+i, b_off, 0, -1); //-B


        circuit.bin_gates.emplace_back(block_len+len+i, sig_off,sig_off, -1,0);//-sigma^2
        circuit.uni_gates.emplace_back(block_len+len+i, sig_off, 0, 1);  //sigma
        circuit.uni_gates.emplace_back(block_len+len+i, b_off, 0, 1); //+B

        circuit.uni_gates.emplace_back(block_len+len*2+i, sum_off, 0, -1);  //-SUM
        */
    }
    
    calcNormalLayer(circuit, layer_id);
    for (i64 i = 0; i < block_len; i++)
    {
        assert(!val[layer_id][i].isNegative());  //only need to check the sparse items of these
        assert(!val[layer_id][block_len+i].isNegative());  
        if(i<block_len)
        {
            if(sparsity_map[i]==false)
                assert(val[layer_id][i].isZero());
            if(sparsity_map[i]==false)
                assert(val[layer_id][i+block_len].isZero());
        }
    }
    for (i64 i = 0; i < len; i++)
    {
        int b_off=orgsize+10+block_len+len+i;
        assert(val[layer_id][2*block_len+len+i].isZero());
        assert(val[layer_id][2*block_len+i].isZero());
    }

    layer_id++;
}

void neuralNetwork::ln_checker_layer3(layer &circuit, i64 &layer_id, int ln_id, int ey,int cy,int real_cn_in,bool* sparsity_map)
{
    i64 block_len = len* channel_in;
    int output_size=len*real_cn_in;
    initLayer(circuit, output_size, layerType::LAYER_NORM_3);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=0;
    int orgsize=ln_aux_start;  // use previous layer's org size
    
    for (i64 i = 0; i < len; i++)
    {
        for (i64 co = 0; co < real_cn_in; ++co) 
        {
            int g=matIdx(i,co,channel_in);
            int p=matIdx(i, co, real_cn_in);
            int d2_off=orgsize+10+block_len+len*4+p;
            int c1=orgsize;
            circuit.uni_gates.emplace_back(p, d2_off, 0, -1);  
            circuit.bin_gates.emplace_back(p, g, g+block_len,1,1);
        }
        
    }
    calcNormalLayer(circuit, layer_id);
    for(int g=0;g<output_size;g++)
    {
        assert(val[layer_id][g].isZero());   // assert will fail without UNDEF NDEBUG on cmake
    }
    layer_id++;
    q_offset=ln_aux_start+10;  // get the rounded result for next computation
} 



void neuralNetwork::gelu_checker_layer1(layer &circuit, i64 &layer_id, int real_cn_out, int ea,int ca,int eb,int cb,int ec,int cc,int ed,int cd,int ex,int cx,int ey,int cy)
{
    int m1=multi_max::max(0,-ec,-eb-ex,-ea-2*ex,ex-ey);
    int m2=multi_max::max(0,-ed,-ex);
    ll C1 = cx*(1ll<<(m1));
    ll C2 = (1ll<<(ea +2*ex+m1))*cx*ca*cx*cx;
    ll C3 = (1ll<<(eb+ex+m1))*cx*cb*cx;
    ll    C4 = cx*cc * (1ll<<(ec +m1));
    ll    C5 = (1ll<<(ey-ex+m1))*cy;
    
    ll    C6 = cd * (1ll<<(ed +m2));
    ll    C7 = cx * (1ll<<(ex +m2));

    i64 block_len = len* channel_out;
    int output_size=6*block_len;
    initLayer(circuit, output_size, layerType::GELU_1);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=block_len*2;
    int orgsize=val[0].size() ;
    val[0].resize(orgsize +10 + block_len*3 + len*real_cn_out*3 );// Const; y; abs; t; d1; d2; d3
    for(int i=orgsize;i<val[0].size();i++)
        val[0][i]=0;
    gelu_aux_start=orgsize; 
    positive_check+=len*real_cn_out*3; //add positive check
    total_relu_in_size += 10+ len*real_cn_out*3 + block_len*3;
    for (i64 g = 0; g < block_len; ++g) 
    {
        if(g%channel_out>=real_cn_out)
            continue;
        ll qx=convert(val[0][g+q_offset]);
        
        int y_off=orgsize+10;
        int abs_off=orgsize+10+block_len;
        int t_off=orgsize+10+block_len*2;
        int d1_off=orgsize+10+block_len*3;  
        int d2_off=orgsize+10+block_len*3+ len*real_cn_out;
        int d3_off=orgsize+10+block_len*3+ len*real_cn_out*2;
        
        ll abs,t;
        if(qx<0)
            val[0][abs_off+g]=abs=-qx;
        else
            val[0][abs_off+g]=abs=qx;
        if(C6>=C7*abs)
        {
            val[0][t_off+g]=t=1;
        }
        else
        {
            val[0][t_off+g]=t=0;
        }
        int gp=g/channel_out* real_cn_out+g%channel_out;
        val[0][d1_off+gp]=abs+1;
        val[0][d2_off+gp]=t+(1-2*t)*(C7*abs-C6);
        assert(!val[0][d1_off+gp].isNegative());
        assert(!val[0][d2_off+gp].isNegative());
        double inner=(double)ca*cx*cx*qx*qx*pow(2,ea+2*ex)-cb*cx*abs*pow(2,eb+ex)+cc*pow(2,ec);
        double middle=(double)qx+abs-abs*t*inner;
        double final=(double)cx*pow(2,ex-1-ey)*middle/cy;
        ll y=round(final);
        val[0][y_off+g]=(ll)y;
        ll term1=(2*y + 1)*C5 + 1 - C1 * qx - C1 * abs + C2 *t*abs*abs*abs-C3 *t* abs*abs+C4 *t* abs;
        ll term2= C1 * qx + C1 * abs - C2 *t*abs*abs*abs+C3 *t* abs*abs-C4 *t* abs-(2*y-1)*C5+1;
        assert(term1>0);
        assert(term2>0);
        val[0][d3_off+gp]= Fr(term1)*Fr(term2);
    }
    val[0][orgsize]=1; 
    for (i64 g = 0; g < block_len; ++g) 
    {
        if(g%channel_out>=real_cn_out)
            continue;
        int gp=g/channel_out* real_cn_out+g%channel_out;
        int abs_off=orgsize+10+block_len+g;
        
        int t_off=orgsize+10+block_len*2+g;

        int d1_off=orgsize+10+block_len*3+gp;  
        int d2_off=orgsize+10+block_len*3+ len*real_cn_out+gp;

        int c1=orgsize;

        int q_off=g+q_offset;

        int abs_mult_t_off=g;
        int q_square_off=g+block_len;
        int dt1_off=g+block_len*2;
        int dt2_off=g+block_len*3;
        int t2_off=g+block_len*4;
        int abs_check_off=g+block_len*5;
        circuit.bin_gates.emplace_back(abs_mult_t_off, abs_off, t_off, 1, 0);  // sc, layer
        
        circuit.bin_gates.emplace_back(q_square_off, q_off, q_off, 1, 0);  // sc, layer

        circuit.uni_gates.emplace_back(dt1_off,abs_off,0,1);  //g,u,lu,sc
        circuit.uni_gates.emplace_back(dt1_off,c1,0,1);  //g,u,lu,sc
        circuit.uni_gates.emplace_back(dt1_off,d1_off,0,-1);  //g,u,lu,sc

        circuit.uni_gates.emplace_back(dt2_off,c1,0,-C6); 
        circuit.uni_gates.emplace_back(dt2_off,t_off,0,2ll*C6+1); 
        circuit.uni_gates.emplace_back(dt2_off,abs_off,0,C7); 
        circuit.bin_gates.emplace_back(dt2_off, abs_off, t_off, -2ll*C7, 0);
        circuit.uni_gates.emplace_back(dt2_off,d2_off,0,-1);  //g,u,lu,sc

        circuit.bin_gates.emplace_back(t2_off, t_off, t_off, 1, 0);  // sc, layer
        circuit.uni_gates.emplace_back(t2_off, t_off, 0,-1);  //g,u,lu,sc

        circuit.bin_gates.emplace_back(abs_check_off, q_off, q_off, 1, 0);  // sc, layer
        circuit.bin_gates.emplace_back(abs_check_off, abs_off, abs_off, -1, 0);  // sc, layer
    }
    calcNormalLayer(circuit, layer_id);
    

    for(int g=0;g<block_len;g++)
    {
        //TODO round for python, work wierd
        assert(!val[layer_id][g].isNegative());   // assert will fail without UNDEF NDEBUG on cmake
        assert(!val[layer_id][g+block_len].isNegative());
        if(g%channel_out>=real_cn_out)
        {
            assert(val[layer_id][g].isZero());   
            assert(val[layer_id][g+block_len].isZero());
        }
        assert(val[layer_id][g+block_len*2].isZero());
        assert(val[layer_id][g+block_len*3].isZero());
        assert(val[layer_id][g+block_len*4].isZero());
        assert(val[layer_id][g+block_len*5].isZero());
    }
    layer_id++;
} 

void neuralNetwork::gelu_checker_layer2(layer &circuit, i64 &layer_id, int real_cn_out,int ea,int ca,int eb,int cb,int ec,int cc,int ed,int cd,int ex,int cx,int ey,int cy)
{
    int m1=multi_max::max(0,-ec,-eb-ex,-ea-2*ex,ex-ey);
    int m2=multi_max::max(0,-ed,-ex);
    ll C1 = cx*(1ll<<(m1));
    ll C2 = (1ll<<(ea +2*ex+m1))*cx*ca*cx*cx;
    ll C3 = (1ll<<(eb+ex+m1))*cx*cb*cx;
    ll    C4 = cx*cc * (1ll<<(ec +m1));
    ll    C5 = (1ll<<(ey-ex+m1))*cy;
    
    ll    C6 = cd * (1ll<<(ed +m2));
    ll    C7 = cx * (1ll<<(ex +m2));
    i64 block_len = len* channel_out;
    int output_size=2*block_len;
    initLayer(circuit, output_size, layerType::GELU_2);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    int orgsize=gelu_aux_start;  // use previous layer's org size

    for (i64 g = 0; g < block_len; ++g) 
    {
        if(g%channel_out>=real_cn_out)
            continue;
        int y_off=orgsize+10+g;
        int abs_off=orgsize+10+block_len+g;
        
        int t_off=orgsize+10+block_len*2+g;
        int c1=orgsize;

        int q_off=g+q_offset;

        int term1_off=g;
        int term2_off=g+block_len;

        int abs_mult_t_off=g;
        int q_square_off=g+block_len;

        circuit.uni_gates.emplace_back(term1_off, c1, 0, C5+1);  // C5+1
        circuit.uni_gates.emplace_back(term1_off, y_off, 0, C5*2);  // 2*C5*y
        circuit.uni_gates.emplace_back(term1_off, abs_off, 0, -C1);  
        circuit.uni_gates.emplace_back(term1_off, q_off, 0, -C1);  
        circuit.uni_gates.emplace_back(term1_off, abs_mult_t_off, layer_id-1, C4);  
        circuit.bin_gates.emplace_back(term1_off, abs_mult_t_off,q_square_off,C2,1);
        circuit.bin_gates.emplace_back(term1_off, q_square_off, t_off,-C3,2);


        circuit.uni_gates.emplace_back(term2_off, c1, 0, C5+1);  // C5+1
        circuit.uni_gates.emplace_back(term2_off, y_off, 0, -C5*2);  // 2*C5*y
        circuit.uni_gates.emplace_back(term2_off, abs_off, 0, C1);  
        circuit.uni_gates.emplace_back(term2_off, q_off, 0, C1);  
        circuit.uni_gates.emplace_back(term2_off, abs_mult_t_off, layer_id-1, -C4);  
        circuit.bin_gates.emplace_back(term2_off, abs_mult_t_off,q_square_off,-C2,1);
        circuit.bin_gates.emplace_back(term2_off, q_square_off, t_off,C3,2);
    }
    calcNormalLayer(circuit, layer_id);
    for(int g=0;g<block_len;g++)
    {
        
        if(g%channel_out>=real_cn_out)
        {
            assert(val[layer_id][g].isZero());   
            assert(val[layer_id][g+block_len].isZero());
        }
        assert(!val[layer_id][g].isNegative());   // assert will fail without UNDEF NDEBUG on cmake
        assert(!val[layer_id][g+block_len].isNegative());
    }
    layer_id++;
} 


void neuralNetwork::gelu_checker_layer3(layer &circuit, i64 &layer_id,int real_cn_out, int ea,int ca,int eb,int cb,int ec,int cc,int ed,int cd,int ex,int cx,int ey,int cy)
{

    i64 block_len = len* channel_out;
    int output_size=len* real_cn_out;
    initLayer(circuit, output_size, layerType::GELU_3);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=0;
    int orgsize=gelu_aux_start;  // use previous layer's org size

    for (i64 g = 0; g < block_len; ++g) 
    {
        if(g%channel_out>=real_cn_out)
            continue;
        int gp=g/channel_out* real_cn_out+g%channel_out;

        int d3_off=orgsize+10+block_len*3+2*real_cn_out*len+gp;
        int c1=orgsize;

        circuit.uni_gates.emplace_back(gp, d3_off, 0, -1);  
        circuit.bin_gates.emplace_back(gp, g, g+block_len,1,1);
    }
    calcNormalLayer(circuit, layer_id);
    for(int g=0;g<output_size;g++)
        assert(val[layer_id][g].isZero());   
    layer_id++;
    q_offset=gelu_aux_start+10;  
} 


// we place the computation of after fcon layer here
void neuralNetwork::roundLayer(layer &circuit, i64 &layer_id, float scale,bool* sparsity_map) 
{
    i64 block_len = len* channel_out;
    int c,m;
    pair<int,int> pm=search(scale);
    m=pm.first;
    c=pm.second;
    float virtual_scale=c*pow(2,m);
    i64 size = block_len; 
    initLayer(circuit, size, layerType::RELU);  //TODO: the output dim of the layer
    circuit.need_phase2=true;
    circuit.zero_start_id=0;
    int orgsize=val[0].size() ;
    val[0].resize(orgsize + 20 + block_len*2);// Const; Q; delta
    
    q_offset=orgsize + 20 ;  //TODO set the input offset of the next matrix
    total_relu_in_size += 20+ block_len*2; //TODO: need to update here, for all aux vars added
    val[0][orgsize]=1; 
    int M=max(-m,0);
    for(i64 g = 0; g < block_len; ++g) 
    {
        int qq=g+orgsize+20;
        double fm;
        double fz;
        ll p=convert(val[layer_id-1][g]);
        ll q=round(p*c*pow(2,m));
        val[0][qq]=q;  // compute non-linear round
        int s=qq+block_len;
        
        val[0][s]=Fr(p*c*(1ll<<(m+M+1))+(1<<M)-q*(1<<(M+1)))*Fr(q*(1<<(M+1))+(1<<M)-c*(1ll<<(m+M+1))*p);
        assert(!val[0][s].isNegative());
    }
    for (i64 g = 0; g < block_len; ++g) 
    {
        int p=g;
        int q=g+orgsize+20;
        int c1=orgsize;
        int s=q+block_len;
        if(sparsity_map)
        {
            if(!sparsity_map[g])
                continue;
        }
        circuit.bin_gates.emplace_back(g, p, p, -(1ll<<(2*m+2*M+2))*c*c , 1); //  , 
        circuit.bin_gates.emplace_back(g, p, q, (1ll<<(m+2*M+3))*c ,2);  //
        circuit.bin_gates.emplace_back(g, q, q, -(1<<(2*M+2)) ,0);
        circuit.uni_gates.emplace_back(g, c1, 0, (1ll<<(2*M)) );  // this public input is one
        circuit.uni_gates.emplace_back(g, s, 0, -1);
    }   
    calcNormalLayer(circuit, layer_id);
    for(int i=0;i<block_len;i++)
    {
        int p=i;
        int q=i+orgsize+20;
        assert(val[layer_id][i].isZero());   
    }
    layer_id++;
}


void neuralNetwork::multi_head_matrix_QK(layer &circuit, i64 &layer_id)
{
    const int HEAD=12;
    const int HSIZE=64;
    int output_size=HEAD*len*(len+1)/2;
    initLayer(circuit, output_size, layerType::MHA_QK);
    circuit.need_phase2 = true;
    for(int head=0;head<HEAD;head++)
    {
        int T=0;
        for(int i=0;i<len;i++)
        for(int j=0;j<=i;j++)
        {
            int targ_gate=head*len*(len+1)/2+T;
            for(int k=0;k<HSIZE;k++)
            {
                // (i, head*64+k)
                // (j, HEAD*HSIZE+head*64+k)
                int col_i=head*64+k;
                int col_j=HEAD*HSIZE+head*64+k;
                int idi=i*channel_out+col_i;
                int idj=j*channel_out+col_j;
                int gate_i=q_offset+idi;
                int gate_j=q_offset+idj;
                circuit.bin_gates.emplace_back(targ_gate, gate_i,gate_j, 1,0 );
            }
            ++T;
        }
    }
    calcNormalLayer(circuit, layer_id);
    layer_id++;
}

void neuralNetwork::compute_e_table()
{
    double St=pow(2,-9),Se=pow(2,-20);
    for(int i=0;i<655360;i++)
    {
        int t=round(exp(-St*i)/Se);
        table[i]=max(t,1);  //TODO: avoid sum_Ei=0, occasionally happens
    }
}

void neuralNetwork::softmax_layer_1(layer &circuit, i64 &layer_id,float SQ,float SK,float Sv, float Sy)
{
    const int HEAD=12;
    const int HSIZE=64;
    int orgsize=val[0].size();
    val[0].resize(orgsize+10+HEAD*2*len+4*HEAD*len*(len+1)/2+HEAD*len*HSIZE+len*channel_in);//sumE,pmax,delta1,delta2,t,E,delta3,Y
    for(int i=orgsize;i<val[0].size();i++)
        val[0][i].clear();
    val[0][orgsize]=1;
    softmax_aux_start=orgsize;
    positive_check+=2*HEAD*len*(len+1)/2+HEAD*len*HSIZE; //add positive check for delta1,delta2,delta3
    exp_check+=HEAD*len*(len+1)/2;  //exp check for (t,E) pair
    total_relu_in_size += 10+2*HEAD*len+4*HEAD*len*(len+1)/2+HEAD*len*HSIZE+len*channel_in;
    int output_size=3*HEAD*len*(len+1)/2+HEAD*len+HEAD*len*HSIZE;
    initLayer(circuit, output_size, layerType::SOFTMAX_1);
    circuit.need_phase2 = true;
    circuit.zero_start_id=2*HEAD*len*(len+1)/2+HEAD*len*HSIZE;
    int e1,c1;
    const float St=pow(2,-9);
    const float Se=pow(2,-16);
    pair<int,int> pm=search(SQ*SK/St);
    e1=pm.first;
    c1=pm.second;
    int eprime=max(-e1-1,0);
    for(int head=0;head<HEAD;head++)
    {
        for(int i=1;i<=len;i++)
        {
            int sum_E_offset=orgsize+10+head*len+i-1;
            int pmax_offset=orgsize+10+HEAD*len+head*len+i-1;
            val[0][sum_E_offset]=0;
            //for(int j=0;j<i;j++)
            //{
            //    int offset=head*len*(len+1)/2+T;
            //    val[0][sum_E_offset]+=val[layer_id-1][offset];
            //}
            ll mx=0;
            for(int j=0;j<i;j++)
            {
                int offset=head*len*(len+1)/2+i*(i-1)/2+j;
                ll S=convert(val[layer_id-1][offset]);
                mx=max(mx,S);
            }
            val[0][pmax_offset]=Fr(mx);
            for(int j=0;j<i;j++) // i is length, j is id
            {
                int offset=head*len*(len+1)/2+i*(i-1)/2+j; //global offset
                int dt1_off=orgsize+10+HEAD*2*len+offset;
                int dt2_off=orgsize+10+HEAD*2*len+HEAD*len*(len+1)/2+offset;
                //int dt3_off=orgsize+10+HEAD*2*len+2*HEAD*len*(len+1)/2+offset;
                int t_off=orgsize+10+HEAD*2*len+2*HEAD*len*(len+1)/2+offset;
                int E_off=orgsize+10+HEAD*2*len+3*HEAD*len*(len+1)/2+offset;
                val[0][dt1_off]=val[0][pmax_offset]-val[layer_id-1][offset]; //pmax-pj
                ll pj_=convert(val[0][dt1_off]);
                
                ll tj=round(c1*pow(2,e1+eprime+1)*pj_/pow(2,eprime+1));
                val[0][t_off]=tj;
                assert(tj>=0 && tj<655360);  //TODO change to 65536
                val[0][E_off]=table[tj];
                val[0][dt2_off]=Fr(c1*(1<<(e1+eprime+1))*pj_+(1<<eprime)-tj*(1<<(eprime+1)))*Fr(-c1*(1<<(e1+eprime+1))*pj_+(1<<eprime)+tj*(1<<(eprime+1)));
                val[0][sum_E_offset]+=table[tj];
                ++T;
            }
        }
    }
    for(int head=0;head<HEAD;head++)
    {
        for(int i=0;i<len;i++)
        {
            for(int j=0;j<HSIZE;j++)
            {
                int out_ij=head*len*HSIZE+i*HSIZE+j;
                for(int k=0;k<=i;k++)
                {
                    int Vkj_offset=q_offset+channel_out*k+2*HEAD*HSIZE+head*HSIZE+j; //Vkj, offset   
                    int E_off=orgsize+10+HEAD*2*len+3*HEAD*len*(len+1)/2+head*len*(len+1)/2+(i*(i+1))/2+k;
                    circuit.bin_gates.emplace_back(out_ij, Vkj_offset,E_off,1,0);
                }
            }
        }
        
        for(int i=1;i<=len;i++)
        {
            int pmax_offset=orgsize+10+HEAD*len+head*len+i-1;
            int sum_E_offset=orgsize+10+head*len+i-1;
            for(int j=0;j<i;j++)
            {
                int offset=head*len*(len+1)/2+(i*(i-1))/2+j; //global offset
                int dt1_off=orgsize+10+HEAD*2*len+offset;
                int check_seg1_offset=2*HEAD*len*(len+1)/2+HEAD*len*HSIZE+offset;
                circuit.uni_gates.emplace_back(check_seg1_offset, dt1_off,0,-1); //-(pm-pi)
                circuit.uni_gates.emplace_back(check_seg1_offset, offset,layer_id-1,-1); //-pi
                circuit.uni_gates.emplace_back(check_seg1_offset, pmax_offset,0,1); //+pmax
                int term1_offset=HEAD*len*HSIZE+offset;
                int term2_offset=HEAD*len*(len+1)/2+HEAD*len*HSIZE+offset;
                int t_off=orgsize+10+HEAD*2*len+2*HEAD*len*(len+1)/2+offset;
                circuit.uni_gates.emplace_back(term1_offset, dt1_off,0,c1*(1<<(e1+1+eprime))); //c1*2^(e1+1+e')
                circuit.uni_gates.emplace_back(term1_offset, orgsize,0,1<<eprime); //+2^e'
                circuit.uni_gates.emplace_back(term1_offset, t_off,0,-(1<<(eprime+1))); //-2^(e'+1)*ti

                circuit.uni_gates.emplace_back(term2_offset, dt1_off,0,-c1*(1<<(e1+1+eprime))); //c1*2^(e1+1+e')
                circuit.uni_gates.emplace_back(term2_offset, orgsize,0,1<<eprime); //+2^e'
                circuit.uni_gates.emplace_back(term2_offset, t_off,0,(1<<(eprime+1))); //-2^(e'+1)*ti
            }
            int sum_e_check=3*HEAD*len*(len+1)/2+HEAD*len*HSIZE+head*len+i-1;
            circuit.uni_gates.emplace_back(sum_e_check, sum_E_offset,0,-1); //-2^(e'+1)*ti
            for(int j=0;j<i;j++) // i is length, j is id
            {
                int offset=head*len*(len+1)/2+(i*(i-1))/2+j; //global offset
                int E_off=orgsize+10+HEAD*2*len+3*HEAD*len*(len+1)/2+offset;
                circuit.uni_gates.emplace_back(sum_e_check, E_off,0,1);
            }
        }
    }
    
    calcNormalLayer(circuit, layer_id);
    for(int i=circuit.zero_start_id;i<val[layer_id].size();i++)
    {
        assert(val[layer_id][i].isZero());
    }
    for(int i=HEAD*len*HSIZE;i<circuit.zero_start_id;i++)
    {
        assert(!val[layer_id][i].isNegative());
    }
    layer_id++;
}

void neuralNetwork::softmax_layer_2(layer &circuit, i64 &layer_id,float SQ,float SK,float Sv, float Sy)
{
    const int HEAD=12;
    const int HSIZE=64;
    int orgsize=softmax_aux_start;
    
    int output_size=2*HEAD*len*HSIZE+HEAD*len*(len+1)/2; //delta3_term1, delta3_term2, delta2_check 
    initLayer(circuit, output_size, layerType::SOFTMAX_2);
    circuit.need_phase2 = true;
    circuit.zero_start_id=2*HEAD*len*HSIZE;
    int e1,c1;
    pair<int,int> pm=search(Sv/Sy);
    e1=pm.first;
    c1=pm.second;
    int eprime=max(-e1,1);
    ll f=0,f2=0;
    for(int head=0;head<HEAD;head++)
    {
        for(int i=0;i<len;i++)
        {
            int sum_E_offset=orgsize+10+head*len+i;
            ll sumE=convert(val[0][sum_E_offset]);
            assert(sumE!=0);
            for(int j=0;j<HSIZE;j++)
            {
                int out_ij=head*len*HSIZE+i*HSIZE+j; //on layer_id-1
                int term1_off=out_ij;
                int term2_off=out_ij+HEAD*len*HSIZE;
                int s_ij=orgsize+10+HEAD*2*len+4*HEAD*len*(len+1)/2+HEAD*len*HSIZE+i*channel_in+head*HSIZE+j;
                int d3_off=orgsize+10+HEAD*2*len+4*HEAD*len*(len+1)/2+head*len*HSIZE+i*HSIZE+j;
                ll Qij=convert(val[layer_id-1][out_ij]);
                ll S=(ll)round(Qij*c1*pow(2,e1)/sumE);
                val[0][s_ij]=S;
                ll d1=(sumE*(1<<(eprime-1))+Qij*c1*(1<<(eprime+e1))-S*sumE*(1<<eprime));
                ll d2=(sumE*(1<<(eprime-1))-Qij*c1*(1<<(eprime+e1))+S*sumE*(1<<eprime));
                val[0][d3_off]=Fr(sumE*(1<<(eprime-1))+Qij*c1*(1<<(eprime+e1))-S*sumE*(1<<eprime))*Fr((sumE*(1<<(eprime-1))-Qij*c1*(1<<(eprime+e1))+S*sumE*(1<<eprime)));
                f2=min(f2,S);
                circuit.uni_gates.emplace_back(term1_off, out_ij,layer_id-1,c1*(1ll<<(e1+eprime)));
                circuit.uni_gates.emplace_back(term1_off, sum_E_offset,0,1ll<<(eprime-1));

                circuit.bin_gates.emplace_back(term1_off, s_ij, sum_E_offset,-(1<<eprime),0);
                circuit.uni_gates.emplace_back(term2_off, out_ij,layer_id-1,-c1*(1ll<<(e1+eprime)));
                circuit.uni_gates.emplace_back(term2_off, sum_E_offset,0,1ll<<(eprime-1));
                circuit.bin_gates.emplace_back(term2_off, s_ij, sum_E_offset,(1<<eprime),0);
            }
        }
    }
    for(int head=0;head<HEAD;head++)
    {        
        for(int i=1;i<=len;i++)
        {
            for(int j=0;j<i;j++)
            {
                int offset=head*len*(len+1)/2+i*(i-1)/2+j; //global offset
                int dt2_off=orgsize+10+HEAD*2*len+HEAD*len*(len+1)/2+offset;
                int term1_offset=HEAD*len*HSIZE+offset;
                int term2_offset=HEAD*len*(len+1)/2+HEAD*len*HSIZE+offset;
                int now_off=2*HEAD*len*HSIZE+offset;
                circuit.bin_gates.emplace_back(now_off, term1_offset, term2_offset,1,1);
                circuit.uni_gates.emplace_back(now_off,dt2_off,0,-1);
            }
        }
    }
    
    
    calcNormalLayer(circuit, layer_id);
    for(int i=circuit.zero_start_id;i<val[layer_id].size();i++)
    {
        assert(val[layer_id][i].isZero());
    }
    for(int i=0;i<circuit.zero_start_id;i++)
    {
        assert(!val[layer_id][i].isNegative());
    }
    layer_id++;
    
}
void neuralNetwork::softmax_layer_3(layer &circuit, i64 &layer_id,float SQ,float SK,float Sv, float Sy)
{
    const int HEAD=12;
    const int HSIZE=64;
    int orgsize=softmax_aux_start;
    
    int output_size=HEAD*len*HSIZE; //delta3_check
    initLayer(circuit, output_size, layerType::SOFTMAX_3);
    circuit.need_phase2 = true;
    circuit.zero_start_id=0;
    for(int head=0;head<HEAD;head++)
    {
        for(int i=0;i<len;i++)
        {
            int sum_E_offset=orgsize+10+head*len+i;
            ll sumE=convert(val[0][sum_E_offset]);
            assert(sumE!=0);
            for(int j=0;j<HSIZE;j++)
            {
                int now=head*len*HSIZE+i*HSIZE+j; //on layer_id-1
                int term1_off=now;
                int term2_off=now+HEAD*len*HSIZE;
                int d3_off=orgsize+10+HEAD*2*len+4*HEAD*len*(len+1)/2+head*len*HSIZE+i*HSIZE+j;
                
                
                circuit.bin_gates.emplace_back(now,term1_off, term2_off, 1,1);
                circuit.uni_gates.emplace_back(now,d3_off ,0,-1);
            }
        }
    }

    calcNormalLayer(circuit, layer_id);
    for(int i=circuit.zero_start_id;i<val[layer_id].size();i++)
    {
        assert(val[layer_id][i].isZero());
    }
    layer_id++;
    q_offset=orgsize+10+HEAD*2*len+4*HEAD*len*(len+1)/2+HEAD*len*HSIZE; 
}

void neuralNetwork::fullyConnLayer(layer &circuit, i64 &layer_id, i64 first_fc_id,  int x_offset, int x_layer) 
{
    i64 size = channel_out*len;
    initLayer(circuit, size, layerType::FCONN);
    circuit.need_phase2 = true;
    val[layer_id].resize(circuit.size);
    for (i64 i = 0; i < len; i++)
    {
        for (i64 co = 0; co < channel_out; ++co) 
        {
            i64 g = matIdx(i, co, channel_out);
            val[layer_id][g]=0;
            //circuit.uni_gates.emplace_back(g, first_bias_id + co, 0, 1);  // our protocol doesn't support adding bias for simplicity
            for (i64 ci = 0; ci < channel_in; ++ci) 
            {
                i64 u = x_offset+matIdx(i, ci, channel_in);
                i64 v = first_fc_id + matIdx(co, ci, channel_in);  // the matrix is distributed as (i,ci)*(co,ci)
                val[layer_id][g]+=val[x_layer][u]*val[0][v];
            }
        }
    }
    layer_id++;
}


void neuralNetwork::refreshFCParam(const fconKernel &fc) {
    channel_in = fc.channel_in;
    channel_out = fc.channel_out;
}

i64 neuralNetwork::getFFTLen() const {
    return 1L << getFFTBitLen();
}

i8 neuralNetwork::getFFTBitLen() const {
    return 0;
}









void neuralNetwork::calcSizeAfterPool(const poolKernel &p) {
}

void neuralNetwork::calcInputLayer(layer &circuit) 
{
    val[0].resize(circuit.size);

    assert(val[0].size() == total_in_size);
    auto val_0 = val[0].begin();

    double num, mx = -10000, mn = 10000;
    vector<double> input_dat;
    int hidden=768;
    for (i64 i=0;i<len;i++)
    {
        for(i64 j=0;j<hidden;j++)
        {
            in >> num; 
            input_dat.push_back(num);
            mx = max(mx, num);
            mn = min(mn, num);
        }
    }
    pair<int,int> pm=search(0.01);  
    input_e=pm.first;
    input_c=pm.second;
    
    double sc=input_c*pow(2,input_e);
    int k=0;
    for (i64 i=0;i<len;i++)
    {
        for(i64 j=0;j<hidden;j++)
        {
            ll s=input_dat[k++]/sc;
            val[0][i*1024+j] = F(s);
        }
        for(i64 j=hidden;j<1024;j++)
            val[0][i*1024+j] =0;
    }

    val_0=val[0].begin()+len*1024;
    for (; val_0 < val[0].begin() + circuit.size; ++val_0) 
        val_0 -> clear();
}



void neuralNetwork::readBias(i64 first_bias_id) {
    auto val_0 = val[0].begin() + first_bias_id;

    double num, mx = -10000, mn = 10000;
    vector<double> input_dat;
    for (i64 co = 0; co < channel_out; ++co) 
    {
        in >> num;
        input_dat.push_back(num);
        mx = max(mx, num);
        mn = min(mn, num);
    }

    for (double i : input_dat)  
        *val_0++ = F((i64) (i * exp2(w_bit + x_bit)));

}

void neuralNetwork::readFconWeight(i64 first_fc_id,int real_r,int real_c,int id) 
{
    double num, mx = -10000, mn = 10000;
    auto val_0 = val[0].begin() + first_fc_id;
    mat_values[id]=new int[4096*1024];
    for (i64 co = 0; co < channel_out; ++co)
        for (i64 ci = 0; ci < channel_in; ++ci) 
        {
            if(co<real_c && ci<real_r)
            {
                mat_values[id][co*channel_in+ci]=rand()%1024;
                val_0[co*channel_in+ci]=mat_values[id][co*channel_in+ci];
            }
            else
            {
                mat_values[id][co*channel_in+ci]=0;
                val_0[co*channel_in+ci]=0;
            }
        }
}

void neuralNetwork::prepareDecmpBit(i64 layer_id, i64 idx, i64 dcmp_id, i64 bit_shift) {
    auto data = abs(val[layer_id].at(idx).getInt64());
    val[0].at(dcmp_id) = (data >> bit_shift) & 1;
}

void neuralNetwork::prepareFieldBit(const F &data, i64 dcmp_id, i64 bit_shift) {
    auto tmp = abs(data.getInt64());
    val[0].at(dcmp_id) = (tmp >> bit_shift) & 1;
}

void neuralNetwork::prepareSignBit(i64 layer_id, i64 idx, i64 dcmp_id) {
    val[0].at(dcmp_id) = val[layer_id].at(idx).isNegative() ? F_ONE : F_ZERO;
}

void neuralNetwork::prepareMax(i64 layer_id, i64 idx, i64 max_id) {
    auto data = val[layer_id].at(idx).isNegative() ? F_ZERO : val[layer_id].at(idx);
    if (data > val[0].at(max_id)) val[0].at(max_id) = data;
}

void neuralNetwork::calcNormalLayer(const layer &circuit, i64 layer_id,bool output) 
{
    val[layer_id].resize(circuit.size);
    for (auto &x: val[layer_id]) 
        x.clear();
    for (auto &gate: circuit.uni_gates) 
    {
        val[layer_id].at(gate.g) = val[layer_id].at(gate.g) + val[gate.lu].at(gate.u) * gate.sc;
    }

    for (auto &gate: circuit.bin_gates) 
    {
        u8 bin_lu = gate.getLayerIdU(layer_id), bin_lv = gate.getLayerIdV(layer_id);
        
        val[layer_id].at(gate.g) = val[layer_id].at(gate.g) + val[bin_lu].at(gate.u) * val[bin_lv][gate.v] * gate.sc;
    }
}

void neuralNetwork::checkNormalLayer(const layer &circuit, i64 layer_id,const vector<vector<F> > & val) 
{
    vector<F> valp;

    valp.resize(val[layer_id].size());
    
    for (auto &x: valp) 
        x.clear();
    for (auto &gate: circuit.uni_gates) 
    {
        assert(gate.g>=0 && gate.g<valp.size());
        assert(gate.u>=0 && gate.u<val[gate.lu].size());
        valp.at(gate.g) += val[gate.lu].at(gate.u) * gate.sc;
    }
    for (auto &gate: circuit.bin_gates) 
    {
        u8 bin_lu = gate.getLayerIdU(layer_id), bin_lv = gate.getLayerIdV(layer_id);
        assert(gate.g>=0 && gate.g<valp.size());
        valp.at(gate.g)+=  val[bin_lu].at(gate.u) * val[bin_lv][gate.v] * gate.sc;
    }
    for(int i=0;i<circuit.size;i++)
        assert(valp[i]==val[layer_id][i]);
}

void neuralNetwork::calcDotProdLayer(const layer &circuit, i64 layer_id) {
    val[layer_id].resize(circuit.size);
    for (int i = 0; i < circuit.size; ++i) val[layer_id][i].clear();

    char fft_bit = circuit.fft_bit_length;
    u32 fft_len = 1 << fft_bit;
    u8 l = layer_id - 1;
    for (auto &gate: circuit.bin_gates)
        for (int s = 0; s < fft_len; ++s)
            val[layer_id][gate.g << fft_bit | s] = val[layer_id][gate.g << fft_bit | s] +
                    val[l][gate.u << fft_bit | s] * val[l][gate.v << fft_bit | s];
}


int neuralNetwork::getNextBit(int layer_id) {
    F mx = F_ZERO, mn = F_ZERO;
    for (const auto &x: val[layer_id]) {
        if (!x.isNegative()) mx = max(mx, x);
        else mn = max(mn, -x);
    }
    i64 x = (mx + mn).getInt64();
    double real_scale = x / exp2(x_bit + w_bit);
    int res = (int) log2( ((1 << (Q - 1)) - 1) / real_scale );
    return res;
}

void neuralNetwork::printLayerValues(prover &pr) {
    for (i64 i = 0; i < SIZE; ++i) 
    {
        for (i64 j = 0; j < std::min(200u, pr.C.circuit[i].size); ++j)
            if (!pr.val[i][j].isZero()) cerr << pr.val[i][j] << ' ';
        cerr << endl;
        for (i64 j = pr.C.circuit[i].zero_start_id; j < pr.C.circuit[i].size; ++j)
            if (pr.val[i].at(j) != F_ZERO) 
            {
                exit(EXIT_FAILURE);
            }
    }
}

void neuralNetwork::printInfer(prover &pr) {
    // output the inference result with the size of (pic_parallel x n_class)
    if (out.is_open()) 
    {
        int n_class = full_conn.back().channel_out;
        for (int p = 0; p < pic_parallel; ++p) {
            int k = -1;
            F v;
            for (int c = 0; c < n_class; ++c) {
                auto tmp = val[SIZE - 1].at(matIdx(p, c, n_class));
                if (!tmp.isNegative() && (k == -1 || v < tmp)) {
                    k = c;
                    v = tmp;
                }
            }
            out << k << endl;
        }
    }
    out.close();
}