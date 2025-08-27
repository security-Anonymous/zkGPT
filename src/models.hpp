#ifndef ZKGPT_HPP
#define ZKGPT_HPP

#include "neuralNetwork.hpp"



class LLM: public neuralNetwork {

public:
    explicit LLM(int depth, int headnum=0, int headdim=0, int attn_dim=0, int linear_dim=0);
};

#endif 
