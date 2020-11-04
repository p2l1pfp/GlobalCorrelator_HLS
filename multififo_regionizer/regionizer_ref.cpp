#include "firmware/regionizer.h"

#include <list>
#include <vector>

template<typename P>
P rePhiObj(const P & t, int phi) {
    P ret = t;
    ret.hwPhi = phi;
    return ret;
}



struct RegionBufferTK { 
    std::list<TkObj> fifos[NTKFIFOS]; 
    TkObj staging_area[NTKFIFOS/2], queue[NTKFIFOS/2];
    void flush() { 
        for (int j = 0; j < NTKFIFOS; ++j) fifos[j].clear(); 
        for (int j = 0; j < NTKFIFOS/2; ++j) clear(staging_area[j]);
        for (int j = 0; j < NTKFIFOS/2; ++j) clear(queue[j]);
    }
    void push(int f, const TkObj & tk, int phi_shift=0) {
        fifos[f].push_front(phiShifted(tk, phi_shift));
    }
    TkObj pop_next() {
        TkObj ret; clear(ret);
        // shift data from each pair of fifos to the staging area
        for (int j = 0; j < NTKFIFOS/2; ++j) {
            if (staging_area[j].hwPt != 0) continue;
            if (!fifos[2*j].empty()) {
                staging_area[j] = fifos[2*j].back();
                fifos[2*j].pop_back(); 
            } else if (!fifos[2*j+1].empty()) {
                staging_area[j] = fifos[2*j+1].back();
                fifos[2*j+1].pop_back(); 
            }
        }
        // then from staging area to output
        for (int j = 0; j < NTKFIFOS/2; ++j) {
            if (staging_area[j].hwPt != 0 && queue[j].hwPt == 0) {
                queue[j] = staging_area[j];
                clear(staging_area[j]);
            }
        }
        for (int j = 0; j < NTKFIFOS/2; ++j) {
            if (queue[j].hwPt != 0) {
                ret = queue[j];
                clear(queue[j]);
                break;
            }
        }
        return ret;
    }
    void pop_all(TkObj out[]) {
        for (int j = 0; j < NTKFIFOS; ++j) {
            if (!fifos[j].empty()) {
                out[j] = fifos[j].back();
                fifos[j].pop_back(); 
            } else {
                clear(out[j]);
            }
        }
    }

};

struct RegionBufferCalo { 
    static const int REGION_SIZE = PFREGION_PHI_SIZE;
    static const int SECTOR_SIZE = REGION_SIZE * 3;
    static const int INT_PI  = (REGION_SIZE * 9)/2;

    RegionBufferCalo() :
        ireg(999), nfifo(0), phi_center(0),
        fifos(), staging_area(), queue(), queue2()
    {
    }    
    
    RegionBufferCalo(unsigned int iregion) :
        ireg(iregion),
        nfifo((iregion % 3 == 0) ? 4 : 8),
        phi_center(iregion * REGION_SIZE),
        fifos(nfifo),
        staging_area(nfifo/2), queue(nfifo/2),  staging_area2(nfifo/4), queue2(nfifo/4)
    {
    }
    unsigned int ireg, nfifo, phi_center;
    std::vector<std::list<HadCaloObj>> fifos; 
    std::vector<HadCaloObj> staging_area, queue, staging_area2, queue2;

    void flush() { 
        for (auto & f : fifos) f.clear(); 
        for (auto & t : staging_area) clear(t);
        for (auto & t : queue) clear(t);
        for (auto & t : staging_area2) clear(t);
        for (auto & t : queue2) clear(t);
    }
    void maybe_push(unsigned int sector, unsigned int fiber, const HadCaloObj & tk) {
        int phi_shift = int(sector) * SECTOR_SIZE - phi_center;
        int local_phi = tk.hwPhi.to_int() + phi_shift;
        if (local_phi >= INT_PI) local_phi -= 2*INT_PI;
        if (local_phi < -INT_PI) local_phi += 2*INT_PI;
        if (std::abs(local_phi) <= REGION_SIZE/2+PFREGION_PHI_BORDER) {
            int ififo = fiber + ((sector == ireg/3) ? 0 : 4);
            //if (fiber == 0) printf("test calo sec %u -> reg %u: phi calo %+4d  global %+4d  local %+4d -> ififo %d\n",
            //                            sector, ireg, tk.hwPhi.to_int(), tk.hwPhi.to_int() + int(sector) * SECTOR_SIZE, local_phi, ififo);
            assert(ififo < nfifo);
            fifos[ififo].push_front(rePhiObj(tk, local_phi)); // don't use phiShifted that that has no wrap-around
        }
        //else if (fiber == 0) printf("test calo sec %u -> reg %u: phi calo %+4d  global %+4d  local %+4d -> not accepted\n",
        //                                sector, ireg, tk.hwPhi.to_int(), tk.hwPhi.to_int() + int(sector) * SECTOR_SIZE, local_phi);
    }
    HadCaloObj pop_next() {
        HadCaloObj ret; clear(ret);
        // shift data from each pair of fifos to the staging area
        for (int j = 0; j < nfifo/2; ++j) {
            if (staging_area[j].hwPt != 0) continue;
            for (int i = 2*j; i <= 2*j+1; ++i) {
                if (!fifos[i].empty()) {
                    staging_area[j] = fifos[i].back();
                    fifos[i].pop_back(); 
                    break;
                }
            }
        }
        // then from staging area to queue
        for (int j = 0; j < nfifo/2; ++j) {
            if (staging_area[j].hwPt != 0 && queue[j].hwPt == 0) {
                queue[j] = staging_area[j];
                clear(staging_area[j]);
            }
        }
        // then from queue to staging2
        for (int j = 0; j < nfifo/4; ++j) {
            if (staging_area2[j].hwPt != 0) continue;
            for (int i = 2*j; i <= 2*j+1; ++i) {
                if (queue[i].hwPt != 0) {
                    staging_area2[j] = queue[i];
                    clear(queue[i]);
                    break;
                }
            }
        }
        // then from staging2 to queue2
        for (int j = 0; j < nfifo/4; ++j) {
            if (staging_area2[j].hwPt != 0 && queue2[j].hwPt == 0) {
                queue2[j] = staging_area2[j];
                clear(staging_area2[j]);
            }
        }
        // and finally out
        for (int j = 0; j < nfifo/4; ++j) {
            if (queue2[j].hwPt != 0) {
                ret = queue2[j];
                clear(queue2[j]);
                break;
            }
        }
        return ret;
    }
    void pop_all(HadCaloObj out[]) {
        for (int j = 0; j < nfifo; ++j) {
            if (!fifos[j].empty()) {
                out[j] = fifos[j].back();
                fifos[j].pop_back(); 
            } else {
                clear(out[j]);
            }
        }
    }

};

struct RegionBufferMu { 
    static const int REGION_PHI_HALFSIZE = PFREGION_PHI_SIZE/2+PFREGION_PHI_BORDER;
    static const int REGION_ETA_HALFSIZE = PFREGION_ETA_SIZE/2+PFREGION_ETA_BORDER;
    static const int INT_PI  = (PFREGION_PHI_SIZE * 9)/2;
    std::list<MuObj> fifos[NMUFIBERS]; 
    int etaCenter, phiCenter;
    void flush() { 
        for (int j = 0; j < NMUFIBERS; ++j) fifos[j].clear(); 
    }
    void maybe_push(unsigned int ifiber, const GlbMuObj & gmu) {
        int local_phi = gmu.hwPhi.to_int() - phiCenter;
        int local_eta = gmu.hwEta.to_int() - etaCenter;
        if (local_phi >= INT_PI) local_phi -= 2*INT_PI;
        if (local_phi < -INT_PI) local_phi += 2*INT_PI;
        //printf("try push mu ipt %4d  glb eta %+4d phi %+4d -> local  eta %+4d phi %+4d \n",
        //             gmu.hwPt.to_int(), gmu.hwEta.to_int(),  gmu.hwPhi.to_int(), local_eta, local_phi);
        fflush(stdout);
        if (std::abs(local_phi) <= REGION_PHI_HALFSIZE &&
            std::abs(local_eta) <= REGION_ETA_HALFSIZE) {
            MuObj lmu; 
            lmu.hwPt    = gmu.hwPt;
            lmu.hwPtErr = gmu.hwPtErr;
            lmu.hwEta   = local_eta;
            lmu.hwPhi   = local_phi;
            fifos[ifiber].push_front(lmu); 
        }
    }
    MuObj pop_next() {
        MuObj ret; clear(ret);
        for (int j = 0; j < NMUFIBERS; ++j) {
            if (!fifos[j].empty()) {
                ret = fifos[j].back();
                fifos[j].pop_back();
                break;
            }
        }
        return ret;
    }
    void pop_all(MuObj out[]) {
        for (int j = 0; j < NMUFIBERS; ++j) {
            if (!fifos[j].empty()) {
                out[j] = fifos[j].back();
                fifos[j].pop_back(); 
            } else {
                clear(out[j]);
            }
        }
    }

};

template<typename T, unsigned int NSORT>
struct RegionBuilder {
    T sortbuffer[NSORT];
    void push(bool newevt, const T & in, T outsorted[NSORT]) {
        if (newevt) {
            for (int i = 0; i < NSORT; ++i) { 
                outsorted[i] = sortbuffer[i]; 
                clear(sortbuffer[i]); 
            }
        }
        int i = 0; T work = in;
        while (i < NSORT && in.hwPt <= sortbuffer[i].hwPt) i++;
        while (i < NSORT) { std::swap(work, sortbuffer[i]); i++; } 
    }
    void dump(bool newline=true) {
            printf("buff %p", &sortbuffer[0]);
            for (int i = 0; i < NSORT; ++i) printf(" %3d.%03d", sortbuffer[i].hwPt.to_int(), sortbuffer[i].hwEta.to_int());
            if (newline) printf("\n");
    }
};

template<typename T, unsigned int NSORT, unsigned int NOUT>
struct RegionMux {
    T buffer[NPFREGIONS][NSORT];
    unsigned int iter, ireg;
    RegionMux() { 
        iter = 0; ireg = NPFREGIONS; 
        for (int i = 0; i < NPFREGIONS; ++i) {
            for (int j = 0; j < NSORT; ++j) clear(buffer[i][j]);
        }
    }
    bool stream(bool newevt, T stream_out[NOUT]) {
        if (newevt) { iter = 0; ireg = 0; }
        if (ireg < NPFREGIONS) {
#if defined(ROUTER_NOSTREAM)
            assert(NOUT == NSORT);
            for (int i = 0; i < NOUT; ++i) {
                stream_out[i] = buffer[ireg][i];
            }
#else
            for (int i = 0; i < NOUT; ++i) {
                if (PFLOWII*i < NSORT) {
                    stream_out[i] = buffer[ireg][PFLOWII*i];
                } else {
                    clear(steam_out[i]);
                }
            }
            for (int i = 1; i < NSORT; ++i) {
                buffer[ireg][i-1] = buffer[ireg][i];
            }
#endif
            if (++iter == PFLOWII) {
                ireg++; iter = 0;
            }
            return true;
        } else {
            for (int i = 0; i < NOUT; ++i) {
                clear(stream_out[i]);
            }
            return false;
        }
    }
};

struct RegionizerTK {
    RegionBufferTK buffers[NPFREGIONS];
    RegionBuilder<TkObj,NTKSORTED> builder[NTKSECTORS];
    RegionMux<TkObj,NTKSORTED,NTKOUT> bigmux;
    unsigned int nevt;
    RegionizerTK() { nevt = 0; }
    void flush() { 
        for (int i = 0; i < NPFREGIONS; ++i) buffers[i].flush();
    }
    void read_in(const TkObj in[NTKSECTORS][NTKFIBERS]) {
        for (int i = 0; i < NTKSECTORS; ++i) {
            for (int j = 0; j < NTKFIBERS; ++j) {
                const TkObj & tk = in[i][j];
                if (tk.hwPt == 0) continue;
                buffers[i].push(j, tk);
                int inext = (i+1), iprev = i+NTKSECTORS-1;
                bool link_next = tk.hwPhi >= +(PFREGION_PHI_SIZE/2-PFREGION_PHI_BORDER);
                bool link_prev = tk.hwPhi <= -(PFREGION_PHI_SIZE/2-PFREGION_PHI_BORDER);
                if (link_next) buffers[inext%NTKSECTORS].push(j+2, tk, -PFREGION_PHI_SIZE);
                if (link_prev) buffers[iprev%NTKSECTORS].push(j+4, tk, +PFREGION_PHI_SIZE);
            }
        }
    }
    void write_out(TkObj out[NPFREGIONS]) {
        for (int i = 0; i < NTKSECTORS; ++i) {
#if defined(ROUTER_NOMERGE)
            buffers[i].pop_all(&out[i*NTKFIFOS]);
#else
            out[i] = buffers[i].pop_next();
#endif
        }
    }
    bool run(bool newevt, const TkObj in[NTKSECTORS][NTKFIBERS], TkObj out[NTKSORTED]) {
        if (newevt) { flush(); nevt++; }
        read_in(in);
#if defined(ROUTER_MUX) or defined(ROUTER_NOSTREAM)
        TkObj routed[NPFREGIONS];
        write_out(routed);
        for (int i = 0; i < NPFREGIONS; ++i) {
            builder[i].push(newevt, routed[i], &bigmux.buffer[i][0]);
        }
        return bigmux.stream(newevt && (nevt > 1), out);
#else
        write_out(out);
        return true;
#endif
    }
};


struct RegionizerCalo {
    RegionBufferCalo buffers[NPFREGIONS];
    RegionBuilder<HadCaloObj,NCALOSORTED> builder[NPFREGIONS];
    RegionMux<HadCaloObj,NCALOSORTED,NCALOOUT> bigmux;
    unsigned int nevt;
    RegionizerCalo() { 
        for (int r = 0; r < NPFREGIONS; ++r) buffers[r] = RegionBufferCalo(r);
        nevt = 0; 
    }
    void flush() { 
        for (int i = 0; i < NPFREGIONS; ++i) buffers[i].flush();
    }
    void read_in(const HadCaloObj in[NCALOSECTORS][NCALOFIBERS]) {
        for (int i = 0; i < NCALOSECTORS; ++i) {
            for (int j = 0; j < NCALOFIBERS; ++j) {
                if (in[i][j].hwPt == 0) continue;
                for (int r = 0; r < NPFREGIONS; ++r) {
                    buffers[r].maybe_push(i,j,in[i][j]);
                }
            }
        }
    }
    void write_out(HadCaloObj out[NCALOOUT]) {
        for (unsigned int i = 0, offs = 0; i < NPFREGIONS; ++i) {
#if defined(ROUTER_NOMERGE)
            buffers[i].pop_all(&out[offs]);
            offs += buffers[i].nfifo;
#else
            out[i] = buffers[i].pop_next();
#endif
        }
    }
    bool run(bool newevt, const HadCaloObj in[NCALOSECTORS][NCALOFIBERS], HadCaloObj out[NCALOOUT]) {
        if (newevt) { flush(); nevt++; }
        read_in(in);
#if defined(ROUTER_MUX) or defined(ROUTER_NOSTREAM)
        HadCaloObj routed[NPFREGIONS];
        write_out(routed);
        for (int i = 0; i < NPFREGIONS; ++i) {
            builder[i].push(newevt, routed[i], &bigmux.buffer[i][0]);
        }
        return bigmux.stream(newevt && (nevt > 1), out);
#else
        write_out(out);
        return true;
#endif
    }
};

struct RegionizerMu {
    RegionBufferMu buffers[NPFREGIONS];
    RegionBuilder<MuObj,NMUSORTED> builder[NPFREGIONS];
    RegionMux<MuObj,NMUSORTED,NMUOUT> bigmux;
    unsigned int nevt;
    RegionizerMu(const glbeta_t etaCenter) { 
        for (int i = 0; i < NPFREGIONS; ++i) {
            buffers[i].phiCenter = i * PFREGION_PHI_SIZE;
            buffers[i].etaCenter = etaCenter;
        }
        nevt = 0; 
    }
    void flush() { 
        for (int i = 0; i < NPFREGIONS; ++i) buffers[i].flush();
    }
    void read_in(const GlbMuObj in[NMUFIBERS]) {
        for (int j = 0; j < NMUFIBERS; ++j) {
            if (in[j].hwPt == 0) continue;
            for (int r = 0; r < NPFREGIONS; ++r) {
                buffers[r].maybe_push(j,in[j]);
            }
        }
    }
    void write_out(MuObj out[NMUOUT]) {
        for (unsigned int i = 0, offs = 0; i < NPFREGIONS; ++i, offs += NMUFIBERS) {
#if defined(ROUTER_NOMERGE)
            buffers[i].pop_all(&out[offs]);
#else
            out[i] = buffers[i].pop_next();
#endif
        }
    }
    bool run(bool newevt, const GlbMuObj in[NMUFIBERS], MuObj out[NMUOUT]) {
        if (newevt) { flush(); nevt++; }
        read_in(in);
#if defined(ROUTER_MUX) or defined(ROUTER_NOSTREAM)
        MuObj routed[NPFREGIONS];
        write_out(routed);
        for (int i = 0; i < NPFREGIONS; ++i) {
            builder[i].push(newevt, routed[i], &bigmux.buffer[i][0]);
        }
        return bigmux.stream(newevt && (nevt > 1), out);
#else
        write_out(out);
        return true;
#endif
    }
};


bool tk_router_ref(bool newevent, const TkObj tracks_in[NTKSECTORS][NTKFIBERS], TkObj tracks_out[NTKOUT]) {
    static RegionizerTK impl;
    return impl.run(newevent, tracks_in, tracks_out);
}


bool calo_router_ref(bool newevent, const HadCaloObj calo_in[NCALOSECTORS][NCALOFIBERS], HadCaloObj calo_out[NCALOOUT]) {
    static RegionizerCalo impl;
    return impl.run(newevent, calo_in, calo_out);
}

bool mu_router_ref(bool newevent, const glbeta_t etaCenter, const GlbMuObj mu_in[NMUFIBERS], MuObj mu_out[NMUOUT]) {
    static RegionizerMu impl(etaCenter);
    return impl.run(newevent, mu_in, mu_out);
}
