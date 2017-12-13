
template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(EmCaloObj emcalo[N], MP7DataWord data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( emcalo[i].hwPtErr, emcalo[i].hwPt );
        data[2*i+1+OFFS] = ( emcalo[i].hwPhi,   emcalo[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(HadCaloObj hadcalo[N], MP7DataWord data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( hadcalo[i].hwEmPt, hadcalo[i].hwPt );
        data[2*i+1+OFFS] = ( hadcalo[i].hwIsEM, hadcalo[i].hwPhi, hadcalo[i].hwEta );
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(TkObj track[N], MP7DataWord data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( track[i].hwPtErr, track[i].hwPt );
        data[2*i+1+OFFS] = ( track[i].hwZ0, track[i].hwPhi, track[i].hwEta );
        data[2*i+1+OFFS][30] = track[i].hwTightQuality;
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(MuObj mu[N], MP7DataWord data[]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        data[2*i+0+OFFS] = ( mu[i].hwPtErr, mu[i].hwPt );
        data[2*i+1+OFFS] = ( mu[i].hwPhi, mu[i].hwEta );
    }
}


template<unsigned int N, unsigned int OFFS> 
inline void mp7_unpack(MP7DataWord data[], EmCaloObj emcalo[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        emcalo[i].hwPt    = data[2*i+0+OFFS](15, 0);
        emcalo[i].hwPtErr = data[2*i+0+OFFS](31,16);
        emcalo[i].hwEta   = data[2*i+1+OFFS](9,  0);
        emcalo[i].hwPhi   = data[2*i+1+OFFS](19,10);
    }
}

template<unsigned int N, unsigned int OFFS> 
inline void mp7_unpack(MP7DataWord data[], HadCaloObj hadcalo[N]) {
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
inline void mp7_unpack(MP7DataWord data[], TkObj track[N]) {
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
inline void mp7_unpack(MP7DataWord data[], MuObj mu[N]) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
        mu[i].hwPt    = data[2*i+0+OFFS](15, 0);
        mu[i].hwPtErr = data[2*i+0+OFFS](31,16);
        mu[i].hwEta   = data[2*i+1+OFFS](9, 0);
        mu[i].hwPhi   = data[2*i+1+OFFS](19,10);
    }
}


