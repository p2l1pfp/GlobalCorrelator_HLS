#ifndef multififo_regionizer_fifos_h
#define multififo_regionizer_fifos_h

template<unsigned int NBITS=64, unsigned int LEN=64, unsigned int NPTRBITS=6>
class rolling_ram_fifo {
    public: 
        typedef ap_uint<NBITS>   word_t;
        typedef ap_uint<NPTRBITS> ptr_t;
        rolling_ram_fifo() ;
        void update(
                bool roll,        
                const word_t & din,
                bool write,
                word_t & dout,   
                bool & valid,
                bool full,    
                bool &roll_out
        );

    private:
        bool roll_delayed, rolled;
        ptr_t wr_ptr;  // where we have last read the data
        ptr_t rd_ptr;  // where we have read the data
        word_t data[LEN];
        word_t cache;
        bool cache_valid;
};


template<unsigned int NBITS, unsigned int LEN, unsigned int NPTRBITS>
rolling_ram_fifo<NBITS,LEN,NPTRBITS>::rolling_ram_fifo() {
    rd_ptr = 2;
    wr_ptr = 1;
    roll_delayed = 0; rolled = 0;
    cache_valid  = 0;
}


template<unsigned int NBITS, unsigned int LEN, unsigned int NPTRBITS>
void rolling_ram_fifo<NBITS,LEN,NPTRBITS>::update(bool roll,     
                                 const word_t & din, bool write,
                                 word_t & dout, bool & valid, bool full, bool &roll_out) 
{
    #pragma HLS DEPENDENCE variable=data inter false
    #pragma HLS inline

    word_t mem_out = data[rd_ptr]; 
    bool mem_out_valid = (rd_ptr < wr_ptr);

    bool full_and_valid_out = full && cache_valid && !roll_delayed && !rolled;
    if (full_and_valid_out) {
        dout = cache; 
        valid = cache_valid;
    } else {
        cache = mem_out;
        cache_valid = mem_out_valid;
        dout = mem_out;
        valid = mem_out_valid;
    }
    roll_out = roll_delayed;
    rolled = roll_out;

    if (roll) {
        rd_ptr = 1;
    } else if (not(full_and_valid_out) && mem_out_valid) {
        rd_ptr++;
    }

    // implement write port
    if (roll) wr_ptr = 1;
    if (write) {
        data[wr_ptr] = din;
        wr_ptr++;
    }

    roll_delayed = roll;
}

template<unsigned int NBITS=64>
class fifo_merge2_simple {
    public:    
        typedef ap_uint<NBITS> word_t;
        fifo_merge2_simple() ;
        void update(
                bool roll,        
                const word_t & din1,
                const word_t & din2,
                bool valid1, 
                bool valid2,
                word_t & dout,   
                bool  & valid,
                bool &full1,    
                bool &full2,    
                bool &roll_out
        );

    private:
        word_t queue_;
        bool  queue_valid_, full2_;
};

template<unsigned int NBITS>
fifo_merge2_simple<NBITS>::fifo_merge2_simple() {
    queue_valid_ = false; full2_ = false;
}


template<unsigned int NBITS>
void fifo_merge2_simple<NBITS>::update(bool roll,     
                              const word_t & din1,
                              const word_t & din2,
                              bool valid1, 
                              bool valid2,
                              word_t & dout,   
                              bool  & valid,
                              bool  & full1,    
                              bool  & full2,    
                              bool  & roll_out) 
{
    #pragma HLS inline
    if (roll) {
        dout = (valid1 ? din1 : din2);
        valid = valid1 || valid2;
        queue_ = din2;
        queue_valid_ = valid1 && valid2;
        full2_ = valid1 && valid2;
    } else {
        bool load2 = (valid1 || queue_valid_) && !full2_;
        dout = (valid1 ? din1 : (queue_valid_ ? queue_ : din2));
        valid = valid1 || valid2 || queue_valid_;
        full2_ = valid1 && (valid2 || queue_valid_);
        if (load2) {
            queue_ = din2;
            queue_valid_ = valid2;
        } else {
            queue_valid_ = valid1 and queue_valid_;
        }
    }
    roll_out = roll;
    full1 = false;
    full2 = full2_;
}

template<unsigned int NBITS=64>
class fifo_merge2_full {
    public:    
        typedef ap_uint<NBITS> word_t;
        fifo_merge2_full () ;
        void update(
                bool roll,        
                const word_t & din1,
                const word_t & din2,
                bool valid1, 
                bool valid2,
                bool full,
                word_t & dout,   
                bool  & valid,
                bool &full1,    
                bool &full2,    
                bool &roll_out
        );

    private:
        word_t q_[2], out_;
        bool  q_valid_[2], full_[2], valid_, roll_;
};


template<unsigned int NBITS>
fifo_merge2_full<NBITS>::fifo_merge2_full() {
    clear(q_[0]); clear(q_[1]); clear(out_);
    q_valid_[0] = false; full_[0] = false;
    q_valid_[1] = false; full_[1] = false;
    valid_ = false; roll_ = false; 
}

template<unsigned int NBITS>
void fifo_merge2_full<NBITS>::update(bool roll,     
                              const word_t & din1,
                              const word_t & din2,
                              bool valid1, 
                              bool valid2,
                              bool full,
                              word_t & dout,   
                              bool  & valid,
                              bool  & full1,    
                              bool  & full2,    
                              bool  & roll_out) 
{
    #pragma HLS inline
    if (roll) {
        out_ = (valid1 ? din1 : din2);
        valid_ = valid1 || valid2;
        q_[0] = din1;
        q_valid_[0] = false;
        q_[1] = din2;
        q_valid_[1] = valid1 && valid2;
        full_[0] = false;
        full_[1] = valid1 && valid2;
    } else if (full && valid_ && !roll_) {
        if (!full_[0]) {
            q_[0] = din1; q_valid_[0] = valid1; full_[0] = valid1;
        }
        if (!full_[1]) {
            q_[1] = din2; q_valid_[1] = valid2; full_[1] = valid2;
        }
    } else {
        bool load2 = (valid1 || q_valid_[0] || q_valid_[1]) && !full_[1];
        out_ = (q_valid_[0] ? q_[0] : (valid1 ? din1 : (q_valid_[1] ? q_[1] : din2)));
        valid_ = valid1 || valid2 || q_valid_[0] || q_valid_[1];
        full_[0] = false;
        full_[1] = (valid1 || q_valid_[0]) && (valid2 || q_valid_[1]);
        if (load2) {
            q_[1] = din2;
            q_valid_[1] = valid2;
        } else {
            q_valid_[1] = (valid1 || q_valid_[0]) && q_valid_[1];
        }
        q_valid_[0] = false; // important: don't move above the previous if(), which uses q_valid_[0]!
    }
    roll_ = roll;
    roll_out = roll_;
    dout = out_;
    valid = valid_;
    full1 = full_[0];
    full2 = full_[1];
}

template<unsigned int NBITS=64>
class fifo_merge3 {
    public:    
        typedef ap_uint<NBITS> word_t;
        fifo_merge3 () ;
        void update(
                bool roll,        
                const word_t & din1,
                const word_t & din2,
                const word_t & din3,
                bool valid1, 
                bool valid2,
                bool valid3,
                word_t & dout,   
                bool  & valid,
                bool &full1,    
                bool &full2,    
                bool &full3,    
                bool &roll_out
        );
    private:
        word_t q2_, q3_;
        bool  q2_valid_, q3_valid_, full2_, full3_, full2old_, full3old_;
};


template<unsigned int NBITS>
fifo_merge3<NBITS>::fifo_merge3() {
    clear(q2_); clear(q3_);
    q2_valid_ = false; full2_ = false;
    q3_valid_ = false; full3_ = false;
}


template<unsigned int NBITS>
void fifo_merge3<NBITS>::update(bool roll,     
                              const word_t & din1,
                              const word_t & din2,
                              const word_t & din3,
                              bool valid1, 
                              bool valid2,
                              bool valid3,
                              word_t & dout,   
                              bool  & valid,
                              bool  & full1,    
                              bool  & full2,    
                              bool  & full3,    
                              bool  & roll_out) 
{
    #pragma HLS inline
    if (roll) {
        dout = (valid1 ? din1 : (valid2 ? din2 : din3));
        valid = valid1 || valid2 || valid3;
        q2_ = din2;
        q3_ = din3;
        q2_valid_ = valid1 && valid2;
        full2_    = valid1 && valid2;
        q3_valid_ = (valid1 || valid2) && valid3;
        full3_    = (valid1 || valid2) && valid3;
    } else {
        bool load2 = (valid1 || q2_valid_)                        && !full2_;
        bool load3 = (valid1 || valid2 || q2_valid_ || q3_valid_) && !full3_;
        bool q3_valid = q3_valid_, q2_valid = q2_valid_; // cache current values
        full2_ = valid1                          && (valid2 || q2_valid_);
        full3_ = (valid1 || valid2 || q2_valid_) && (valid3 || q3_valid_);
        dout = valid1 ? din1 : (
                 q2_valid ? q2_ : (
                    valid2 ? din2 : (
                        q3_valid ? q3_ :
                            din3)));
        valid = valid1 || valid2 || valid3 || q2_valid || q3_valid_;
        if (load2) {
            q2_ = din2;
            q2_valid_ = valid2;
        } else {
            q2_valid_ = valid1 && q2_valid;
        }
        if (load3) {
            q3_ = din3;
            q3_valid_ = valid3;
        } else {
            q3_valid_ = (valid1 || valid2 || q2_valid) && q3_valid;
        }
    }
    roll_out = roll;
    full1 = false;
    full2 = full2_;
    full3 = full3_;
}

#endif
