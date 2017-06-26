#ifndef RANDOM_INPUTS_H
#define RANDOM_INPUTS_H

class RandomPFInputs {
public:
	RandomPFInputs(int seed) { srand(seed); }
	~RandomPFInputs() { }
	bool nextRegion(CaloObj calo[NCALO], TkObj track[NTRACK], MuObj mu[NMU], z0_t & hwZPV) {
		const float PT_SCALE = 4.0;
		const float ETAPHI_SCALE = (4*180/M_PI);
		const float Z0_SCALE = 20.;

		int ncharged = (rand() % NTRACK/2) + NTRACK/2;
		int nneutral = (rand() % ((3*NCALO)/4));
		for (int i = 1; i < nneutral && i < NCALO; i += 2) {
			float pt = (rand()/float(RAND_MAX))*80+1, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
			calo[i].hwPt  = pt * PT_SCALE;
			calo[i].hwEta = eta * ETAPHI_SCALE;
			calo[i].hwPhi = phi * ETAPHI_SCALE;
		}
		float zPV = (rand()/float(RAND_MAX))*20-10;
		hwZPV = zPV * Z0_SCALE;

		for (int i = 0; i < ncharged && i < NTRACK; ++i) {
			float pt = (rand()/float(RAND_MAX))*50+2, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
			float z = (i % 2 == 0) ? (zPV + (rand()/float(RAND_MAX))*0.7-.35) : ((rand()/float(RAND_MAX))*30-15);
			track[i].hwPt    = pt * PT_SCALE;
			track[i].hwPtErr = (0.2*pt+4) * PT_SCALE; 
			track[i].hwEta = eta * ETAPHI_SCALE;
			track[i].hwPhi = phi * ETAPHI_SCALE;
			track[i].hwZ0  = z * Z0_SCALE;
			int icalo = rand() % NCALO;
			if (i % 3 == 1 || icalo >= NCALO) continue;
			float dpt_calo = ((rand()/float(RAND_MAX))*3-1.5) * (0.2*pt+4);
			float deta_calo = ((rand()/float(RAND_MAX))*0.3-0.15), dphi_calo = ((rand()/float(RAND_MAX))*0.3-0.15);
			if (pt + dpt_calo > 0) {
				calo[icalo].hwPt  += (pt + dpt_calo) * PT_SCALE;
				calo[icalo].hwEta = (eta + deta_calo) * ETAPHI_SCALE;
				calo[icalo].hwPhi = (phi + dphi_calo) * ETAPHI_SCALE;
			}
		}

		for (int i = 0; i < NMU; ++i){
			float pt = (rand()/float(RAND_MAX))*50+2, eta = (rand()/float(RAND_MAX))*2.0-1.0, phi = (rand()/float(RAND_MAX))*2.0-1.0;
			mu[i].hwPt  = pt * PT_SCALE;
			mu[i].hwPtErr = (0.2*pt+4) * PT_SCALE; 
			mu[i].hwEta = eta * ETAPHI_SCALE;
			mu[i].hwPhi = phi * ETAPHI_SCALE;			
		}
		
		return true;
	}
};

#endif
