#ifndef multififo_regionizer_tmux18_utils_h
#define multififo_regionizer_tmux18_utils_h

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>

#include <ap_int.h>


template<typename T, unsigned int TM6CLOCKS>
struct TM18LinkTriplet {
    std::queue<std::pair<T,bool>> links[3];
    
    TM18LinkTriplet() { 
        for (int i = 0; i <  TM6CLOCKS; ++i)  links[1].emplace(T(0), false);
        for (int i = 0; i < 2*TM6CLOCKS; ++i) links[2].emplace(T(0), false);
    };


    template<typename C>
    void push_event(unsigned int iev, const C & objs) {
        auto & q = links[iev % 3];
        unsigned int n = std::min<unsigned int>(3*TM6CLOCKS-3, objs.size()); // let's leave 3 empty frames at the end, so they get reassebled as 1 row of nulls
        //printf("Writing %u objects on link %u for iev %u\n", n, iev%3, iev);
        for (unsigned int i = 0; i < n; ++i) {
            q.emplace(objs[i], true);
        }

        for (unsigned int i = n; i < 3*TM6CLOCKS; ++i) {
            q.emplace(T(0), i < 3*TM6CLOCKS-3);
        }
    }
    template<typename TC, typename BC>
    void pop_frame(TC & values, BC & valids, unsigned int start=0, unsigned int stride=1) {
        for (int i = 0; i < 3; ++i) {
            if (links[i].empty()) { 
                printf("ERROR: link %d is empty (start = %u, stride = %u)\n", i, start, stride); fflush(stdout); 
                continue;
            }
            assert(!links[i].empty());
            auto obj = links[i].front();
            values[start+i*stride] = obj.first;
            valids[start+i*stride] = obj.second;
            links[i].pop();
        }
    }
};

template<typename T, unsigned int TM6CLOCKS>
class TM18LinkMultiplet {

    public:
        TM18LinkMultiplet(unsigned int N) :
            nlinks_(N), links_(N) {}

        template<typename C, typename E>
        void push_links(unsigned int iev, const C objs[], const E & enc) {
            for (unsigned int i = 0; i < nlinks_; ++i) {
                links_[i].push_event(iev, enc(objs[i]));
            }
        }

        template<typename C, typename E>
        void push_link(unsigned int iev, const C objs, const E & enc) {
            assert(nlinks_ == 1);
            links_.front().push_event(iev, enc(objs));
        }


        template<typename TC, typename BC>
        void pop_frame(TC & values, BC & valids, unsigned int start=0, bool group_by_link=true) {
            unsigned int stride = group_by_link ? 1 : nlinks_;
            for (unsigned int i = 0; i < nlinks_; ++i) {
                links_[i].pop_frame(values, valids, start + (group_by_link ? 3*i : i), stride);
            }
        }
    private: 
        unsigned int nlinks_;
        std::vector<TM18LinkTriplet<T,TM6CLOCKS>> links_;
};

class DelayQueue {
    public:
        DelayQueue(unsigned int n) : n_(n), data_(n, 0), ptr_(0) {}
        ap_uint<65> operator()(const ap_uint<65> & in) ;
        void operator()(const ap_uint<64> & in,  const bool & in_valid,
                              ap_uint<64> & out,       bool & out_valid) ;
    private:
        unsigned int n_, ptr_;
        std::vector<ap_uint<65>> data_;
};



#endif
