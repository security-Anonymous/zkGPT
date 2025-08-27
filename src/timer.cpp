#include "timer.hpp"
#include <iostream>
#include <cstring>
using namespace std;

void timer::start() 
{
    assert(status == false);
    t0 = std::chrono::high_resolution_clock::now();
    status = true;
}

void timer::stop(const char* output,bool total,bool restart) 
{
    assert(status == true);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto time_span_sec = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    if(accumulate)
        total_time_sec += time_span_sec.count();
    else
        total_time_sec = time_span_sec.count();
    status = false;
    if(strlen(output)>0)
    {
        if(total)
            cout<<output<<" time: "<<total_time_sec<<endl;
        else
            cout<<output<<" time: "<<time_span_sec.count()<<endl;
    }
    if(restart)
        start();
}
void timer::stop()
{
    stop("",true,false);
} 