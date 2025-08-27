//
// Created by 69029 on 4/12/2021.
//

#undef NDEBUG
#include "circuit.h"
#include "neuralNetwork.hpp"
#include "verifier.hpp"
#include "models.hpp"
#include "global_var.hpp"
#include <iostream>

#include "range_prover.hpp"
#include "hyrax_rp.hpp"
using namespace mcl::bn;
using namespace std;




int main(int argc, char **argv) 
{
    initPairing(mcl::BN254);
    
    // range prover
    range_prover range_prover(12, 12, 64, 768, 2304, 30, 32, 1); // 12 layer, 12 head, 64 channel, 768 head dim, 2304 linear dim, 30 seq len, 32 threads
    range_prover.init();
    range_prover.build();
    double range_prover_time = range_prover.prove();

    // gkr
    prover p;
    LLM nn(12, 12, 64, 768, 2304);  // 12 layer, 12 head, 64 channel, 768 head dim, 2304 linear dim
    nn.create(p, 1);
    verifier v(&p, p.C);
    v.range_prove(range_prover_time);
    v.prove(32); // prove with 32 threads

}

