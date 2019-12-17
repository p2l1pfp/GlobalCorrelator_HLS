template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(EmCaloObj emcalo[N], ap_uint<32> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( emcalo[i].hwPtErr, emcalo[i].hwPt );
        data[2*i+1+OFFS] = ( emcalo[i].hwPhi,   emcalo[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(HadCaloObj hadcalo[N], ap_uint<32> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( hadcalo[i].hwEmPt, hadcalo[i].hwPt );
        data[2*i+1+OFFS] = ( hadcalo[i].hwIsEM, hadcalo[i].hwPhi, hadcalo[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(TkObj track[N], ap_uint<32> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( track[i].hwPtErr, track[i].hwPt );
        data[2*i+1+OFFS] = ( track[i].hwZ0, track[i].hwPhi, track[i].hwEta );
        data[2*i+1+OFFS][30] = track[i].hwTightQuality;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(MuObj mu[N], ap_uint<32> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( mu[i].hwPtErr, mu[i].hwPt );
        data[2*i+1+OFFS] = ( mu[i].hwPhi, mu[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS>
inline void l1pf_pattern_pack(PFChargedObj pfch[N], ap_uint<32> data[]){
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( pfch[i].hwId, pfch[i].hwPt );
        data[2*i+1+OFFS] = ( pfch[i].hwZ0, pfch[i].hwPhi, pfch[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(PFNeutralObj pfne[N], ap_uint<32> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( pfne[i].hwId, pfne[i].hwPt );
        data[2*i+1+OFFS] = ( pfne[i].hwPhi, pfne[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<32> data[], EmCaloObj emcalo[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        emcalo[i].hwPt    = data[2*i+0+OFFS](15, 0);
        emcalo[i].hwPtErr = data[2*i+0+OFFS](31,16);
        emcalo[i].hwEta   = data[2*i+1+OFFS](9,  0);
        emcalo[i].hwPhi   = data[2*i+1+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<32> data[], HadCaloObj hadcalo[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        hadcalo[i].hwPt   = data[2*i+0+OFFS](15, 0);
        hadcalo[i].hwEmPt = data[2*i+0+OFFS](31,16);
        hadcalo[i].hwEta  = data[2*i+1+OFFS](9, 0);
        hadcalo[i].hwPhi  = data[2*i+1+OFFS](19,10);
        hadcalo[i].hwIsEM = data[2*i+1+OFFS][20];
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<32> data[], TkObj track[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        track[i].hwPt    = data[2*i+0+OFFS](15, 0);
        track[i].hwPtErr = data[2*i+0+OFFS](31,16);
        track[i].hwEta   = data[2*i+1+OFFS](9, 0);
        track[i].hwPhi   = data[2*i+1+OFFS](19,10);
        track[i].hwZ0    = data[2*i+1+OFFS](29,20);
        track[i].hwTightQuality = data[2*i+1+OFFS][30];
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<32> data[], MuObj mu[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        mu[i].hwPt    = data[2*i+0+OFFS](15, 0);
        mu[i].hwPtErr = data[2*i+0+OFFS](31,16);
        mu[i].hwEta   = data[2*i+1+OFFS](9, 0);
        mu[i].hwPhi   = data[2*i+1+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS>
inline void l1pf_pattern_unpack(ap_uint<32> data[], PFChargedObj pfch[N]){
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        pfch[i].hwId = data[2*i+0+OFFS](18, 16);
        pfch[i].hwPt = data[2*i+0+OFFS](15, 0);
        pfch[i].hwEta = data[2*i+1+OFFS](9, 0);
        pfch[i].hwPhi = data[2*i+1+OFFS](19,10);
        pfch[i].hwZ0 = data[2*i+1+OFFS](29,20);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<32> data[], PFNeutralObj pfne[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        pfne[i].hwId = data[2*i+0+OFFS](18, 16);
        pfne[i].hwPt = data[2*i+0+OFFS](15, 0);
        pfne[i].hwEta = data[2*i+1+OFFS](9, 0);
        pfne[i].hwPhi = data[2*i+1+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(EmCaloObj emcalo[N], ap_uint<64> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS](31+32, 16+32) = emcalo[i].hwPtErr;
        data[i+OFFS](15+32, 0+32) = emcalo[i].hwPt;
        data[i+OFFS](9, 0) = emcalo[i].hwEta;
        data[i+OFFS](19, 10) = emcalo[i].hwPhi;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(HadCaloObj hadcalo[N], ap_uint<64> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS] = 0;
        data[i+OFFS](31+32,16+32) = hadcalo[i].hwEmPt;
        data[i+OFFS](15+32, 0+32) = hadcalo[i].hwPt;
        data[i+OFFS](9, 0) = hadcalo[i].hwEta;
        data[i+OFFS](19,10) = hadcalo[i].hwPhi;
        data[i+OFFS][20] = hadcalo[i].hwIsEM;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(TkObj track[N], ap_uint<64> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS] = 0;
        data[i+OFFS](31+32,16+32) = track[i].hwPtErr;
        data[i+OFFS](15+32, 0+32) = track[i].hwPt;
        data[i+OFFS](9, 0) = track[i].hwEta;
        data[i+OFFS](19,10) = track[i].hwPhi;
        data[i+OFFS](29,20) = track[i].hwZ0;
        data[i+OFFS][30] = track[i].hwTightQuality;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(MuObj mu[N], ap_uint<64> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS] = 0;
        data[i+OFFS](31+32,16+32) = mu[i].hwPtErr;
        data[i+OFFS](15+32, 0+32) = mu[i].hwPt;
        data[i+OFFS](9, 0) = mu[i].hwEta;
        data[i+OFFS](19,10) = mu[i].hwPhi;
    }
}

template<unsigned int N, unsigned int OFFS>
inline void l1pf_pattern_pack(PFChargedObj pfch[N], ap_uint<64> data[]){
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[i+OFFS] = 0;
        data[i+OFFS](18+32, 16+32) = pfch[i].hwId;
        data[i+OFFS](15+32, 0+32) = pfch[i].hwPt;
        data[i+OFFS](9, 0)  = pfch[i].hwEta;
        data[i+OFFS](19,10) = pfch[i].hwPhi;
        data[i+OFFS](29,20) = pfch[i].hwZ0;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_pack(PFNeutralObj pfne[N], ap_uint<64> data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        //data[i+OFFS](18+32, 16+32) = pfne[i].hwId;
        data[i+OFFS] = 0;
        data[i+OFFS](31+32, 16+32) = pfne[i].hwPtPuppi;
        data[i+OFFS](15+32, 0+32) = pfne[i].hwPt;
        data[i+OFFS](9, 0) = pfne[i].hwEta;
        data[i+OFFS](19,10) = pfne[i].hwPhi;
        data[i+OFFS](22, 20) = pfne[i].hwId;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<64> data[], EmCaloObj emcalo[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        emcalo[i].hwPt    = data[i+OFFS](15+32, 0+32);
        emcalo[i].hwPtErr = data[i+OFFS](31+32,16+32);
        emcalo[i].hwEta   = data[i+OFFS](9,  0);
        emcalo[i].hwPhi   = data[i+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<64> data[], HadCaloObj hadcalo[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        hadcalo[i].hwPt   = data[i+OFFS](15+32, 0+32);
        hadcalo[i].hwEmPt = data[i+OFFS](31+32,16+32);
        hadcalo[i].hwEta  = data[i+OFFS](9, 0);
        hadcalo[i].hwPhi  = data[i+OFFS](19,10);
        hadcalo[i].hwIsEM = data[i+OFFS][20];
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<64> data[], TkObj track[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        track[i].hwPtErr = data[i+OFFS](31+32,16+32);
        track[i].hwPt    = data[i+OFFS](15+32, 0+32);
        track[i].hwEta   = data[i+OFFS](9, 0);
        track[i].hwPhi   = data[i+OFFS](19,10);
        track[i].hwZ0    = data[i+OFFS](29,20);
        track[i].hwTightQuality = data[i+OFFS][30];
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<64> data[], MuObj mu[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        mu[i].hwPt    = data[i+OFFS](15+32, 0+32);
        mu[i].hwPtErr = data[i+OFFS](31+32,16+32);
        mu[i].hwEta   = data[i+OFFS](9, 0);
        mu[i].hwPhi   = data[i+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS>
inline void l1pf_pattern_unpack(ap_uint<64> data[], PFChargedObj pfch[N]){
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        pfch[i].hwId = data[i+OFFS](18+32, 16+32);
        pfch[i].hwPt = data[i+OFFS](15+32, 0+32);
        pfch[i].hwEta = data[i+OFFS](9, 0);
        pfch[i].hwPhi = data[i+OFFS](19,10);
        pfch[i].hwZ0 = data[i+OFFS](29,20);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void l1pf_pattern_unpack(ap_uint<64> data[], PFNeutralObj pfne[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        pfne[i].hwPtPuppi = data[i+OFFS](31+32, 16+32);
        pfne[i].hwPt = data[i+OFFS](15+32, 0+32);
        pfne[i].hwId = data[i+OFFS](22, 20);
        pfne[i].hwEta = data[i+OFFS](9, 0);
        pfne[i].hwPhi = data[i+OFFS](19,10);
    }
}
