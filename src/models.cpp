//
// Created by 69029 on 3/16/2021.
//

#include <tuple>
#include <iostream>
#include "models.hpp"
#include "utils.hpp"
#undef USE_VIRGO




LLM::LLM(int depth, int headnum, int headdim, int attn_dim, int linear_dim):
    neuralNetwork(0,0,0,0, "./data/vgg11/vgg11.cifar.relu-1-images-weights-qint8.csv", "", "", true, headnum, headdim, attn_dim, linear_dim)
{
    len=30;
    pic_parallel=1;
    conv_section.clear();
    positive_check=0;
    exp_check=0;
    layer_num=depth;
    for(int i=0;i<depth;i++)
    {
        full_conn.emplace_back(linear_dim,attn_dim);  //now debugging which act made trouble
        full_conn.emplace_back(attn_dim,attn_dim);
        full_conn.emplace_back(linear_dim,attn_dim);
        full_conn.emplace_back(attn_dim,linear_dim);
    }
    
}