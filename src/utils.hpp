//
// Created by 69029 on 3/9/2021.
//

#ifndef ZKCNN_UTILS_HPP
#define ZKCNN_UTILS_HPP

#include <circuit.h>

int ceilPow2BitLengthSigned(double n);
int floorPow2BitLengthSigned(double n);

char ceilPow2BitLength(u32 n);
char floorPow2BitLength(u32 n);
void field(const char* f,Fr x);



void initBetaTable(vector<F> &beta_g, u8 gLength, const vector<F>::const_iterator &r_0, const vector<F>::const_iterator &r_1,
              const F &alpha, const F &beta,int thd=32);





void initBetaTable(vector<F> &beta_g, u8 gLength, const vector<F>::const_iterator &r, const F &init,int thread=32);

bool check(long x, long y, long nx, long ny);

long matIdx(long x, long y, long n);

long cubIdx(long x, long y, long z, long n, long m);

long tesIdx(long w, long x, long y, long z, long n, long m, long l);

void initLayer(layer &circuit, long size, layerType ty);

long sqr(long x);

double byte2KB(size_t x);

double byte2MB(size_t x);

double byte2GB(size_t x);

F getRootOfUnit(int n);

#endif //ZKCNN_UTILS_HPP
