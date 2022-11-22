

G17B::G17B() {
	memset(this, 0, sizeof(*this));
	gang27 = gang[0];
	aigstop = 0;
	// Init SGUAs tables load permanent data
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int ircol = 0; ircol < 3; ircol++) {//first column
			int* p = tpcol[ircol], dcol = 3 * ist;
			for (int irdig = 0; irdig < 3; irdig++) {//digit1
				for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
					int i = 27 * ist + 9 * ircol + 3 * irdig + irdig2;
					SGUA2& w = tsgua2[i];
					w.i_81 = i; w.i9= dcol + ircol;
					w.col1 = dcol + p[0];		w.col2 = dcol + p[1];
					w.id1 = irdig + 3 * w.col1;
					w.id2 = irdig2 + 3 * w.col2;
				}
			}
		}
	}
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int irdig1 = 0; irdig1 < 3; irdig1++) {//digit1
			for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
				for (int irdig3 = 0; irdig3 < 3; irdig3++) {//digit 2
					int i = 27 * ist + 9 * irdig1 + 3 * irdig2 + irdig3;
					SGUA3& w = tsgua3[i];
					w.i_81 = i;
					w.stack = ist;
					w.col1 = 3 * ist;// minirow first column in gangster
					// pointers to gangster digits
					w.id1 = 9 * ist+irdig1;
					w.id2 = 9 * ist + 3+irdig2;
					w.id3 = 9 * ist + 6 + irdig3;
					if (0) {
						cout << "sgua3 " << w.i_81 << " col1=" << w.col1 + 1
							<< " id1;id2,id3 " << w.id1
							<< ";" << w.id2 << ";" << w.id3 << endl;
					}
				}
			}
		}
	}
}
void G17B::StartInit() {
	memcpy(grid0, genb12.grid0, sizeof grid0);// use first b3
	memcpy(&grid0[54], genb12.bands3[0].band0, sizeof genb12.bands3[0].band0);
	zh2b_g.Init_g0(grid0);// init brute force bands 1+2
	if (op.p1)mincluesb3 = 7;
	else mincluesb3 = 6;
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int* d = gang[i], c = genb12.gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

	//finish sguas2/3 set up with the gangster 
	gsock2.SetAll_0();  gsock3.SetAll_0();
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (register int ist = 0,i=0,stack=07007007; ist < 3; ist++,stack<<=3) {
		//cout << Char27out(stack) << endl;
		for (register int irst = 0; irst < 27; irst++, i++) {//stack
			{//________________________ guas2
				SGUA2& w = tsgua2[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.digs = (1 << w.dig1) | (1 << w.dig2);
				w.valid = 0;
				//cout << i << " " << w.dig1+1 <<" "<< w.dig2+1 << endl;
				//continue;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1] | b3.dpat[w.dig2],
						Bfc = (Bf | (Bf >> 9) | (Bf >> 18)) & 0777;// colummn pattern
					Bf &= stack;
					int nr = 0, bfr = 0, irow;
					if (Bf & 0777) { nr++; bfr |= Bf & 0777; irow = 0; }// row1
					if (Bf & 0777000) { nr++; bfr |= Bf & 0777000; irow = 1; }
					if (Bf & 0777000000) { nr++; bfr |= Bf & 0777000000; irow = 2; }
					if (nr == 1) {// we have a gua2
						//cout << ib3 << " " << i << endl;
						gsock2.setBit(i);
						w.valid = 1;
						b3.g.gsocket2.setBit(i);
						b3.g.pat2[i] = Bf;
						int imini = 3 * irow + ist, mask = 7 << (3 * imini);
						b3.g.ua2_imini[i] = imini;
						int bit27= mask ^ Bf,i27;
						bitscanforward(i27, bit27);
						b3.g.ua2bit27[i] = bit27;
						b3.i_27_to_81[i27] = i;
						b3.i_81_to_27[i] = i27;
					}
					else {// gua46 except if 6 columns
						int ncol = _popcnt32(Bfc);
						if (ncol == 5) {// gua2_4
							b3.g.pat2[i] = bfr;
						}
						else b3.g.pat2[i] = Bf;
					}
				}
			}
			
			
			{ //________________________ guas3

				SGUA3& w = tsgua3[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.dig3 = gang27[w.id3];
				w.digs = 1 << w.dig1 | 1 << w.dig2 | 1 << w.dig3;
				w.valid = 0;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1]
						| b3.dpat[w.dig2] | b3.dpat[w.dig3];// colummn pattern
					Bf &= stack;
					int nr = 0,  irow=0;
					if (Bf & 0777) { nr++;  irow =0; }// row1
					if (Bf & 0777000) { nr++;  irow = 1; }
					if (Bf & 0777000000) { nr++;  irow = 2; }
					if (nr == 1) {// we have a gua3
						gsock3.setBit(i);
						w.valid = 1;
						b3.g.gsocket3.setBit(i);
						b3.g.pat3[i] = Bf;
						int imini = 3 * irow + ist;
						b3.g.ua3_imini[i] = imini;
						b3.i_9_to_81[imini] = i;

					}
				}
			}
		}
	}
}


void G17B::Start() {// processing an entry 
	if (aigstop)return;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;
	StartInit();// do all preliminary setups
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
}

void G17B::StartKnown() {// processing an entry  with a known o follow3
	if (aigstop)return;
	p_cpt2g[0] ++;
	StartInit();// do all preliminary setups
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
}

struct EXTUAS {// get a sorted subset of uas for next "collect" step
	uint64_t  t2[300];
	uint32_t  nt2;
	void GetSortB(uint64_t* t, uint32_t n, uint64_t filter, uint32_t nold = 0) {
		BF128 vsize[25][4];
		uint64_t  t2a[512];
		uint32_t  nt2a = 0;
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		memset(vsize, 0, sizeof vsize);
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			uint64_t cc = w >> 59;
			w &= BIT_SET_2X;
			if (!(w & Fn)) {
				int bloc = nt2a >> 7, ir = nt2a - (bloc << 7);
				vsize[cc][bloc].setBit(ir);
				t2a[nt2a++] = w;// clean the count 
			}
			if (nt2a >= 512)break;
		}
		uint32_t nbl64 = (nt2a + 63) >> 6, x;
		for (int i1 = 4; i1 < 25; i1++) {
			uint64_t* tb64 = vsize[i1]->bf.u64;
			for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
				register uint64_t V = tb64[i2];
				while (bitscanforward64(x, V)) {// switch to 54 mode here
					V ^= (uint64_t)1 << x;
					t2[nt2++] = t2a[x + (i2 << 6)];
					if (nt2 >= 100)return;

				}
			}
		}
	}
	void GetSortA(uint64_t filter, uint32_t nold = 0) {
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			if (!(w & Fn))t2[nt2++] = w;
			if (nt2 >= 100)break;
		}
		if (nt2) {
			if (nt2 > 1)
				for (uint32_t i = 0; i < nt2 - 1; i++) {
					register uint64_t Ri = t2[i];
					for (uint32_t j = i + 1; j < nt2; j++) {
						register uint64_t  Rj = t2[j];
						if ((Ri >> 59) > (Rj >> 59)) { // smaller size first 
							t2[i] = Rj; t2[j] = Ri; Ri = Rj;
						}
					}
				}
			for (uint32_t i = 0; i < nt2; i++)
				t2[i] &= BIT_SET_2X;// clean the count

		}
	}
}extuas;

struct GUAH {// handler for initial collection of guas 2 3
	struct GUA {
		uint64_t tua[128];
		uint32_t nua, type, i81;
		inline void Add(uint64_t ua) {
			if (nua < 128)tua[nua++] = ua;
		}
		inline void AddIf(uint64_t ua) {
			for (uint32_t i = 0; i < nua; i++) {
				if (ua == tua[i]) return;
			}
			tua[nua++] = ua;
		}
		inline void Init(uint32_t n, uint32_t t, uint32_t i) {
			nua = n; type = t; i81 = i;
		}
		int Load(uint64_t* tu, uint64_t bf) {
			int n = 0;
			register uint64_t nbf = ~bf;
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t U = tua[i], cc = _popcnt64(U);
				if (n > 10) break;// 
				if (n > 5 && cc > 12) break;// 
				if (!(U & nbf))tu[n++] = U;
			}
			return n;
		}
		void SortClean() {
			if (nua < 2)return;
			GUA w = *this;
			BF128 vsize[30];
			memset(vsize, 0, sizeof vsize);
			uint64_t tcopy[128];
			memcpy(tcopy, w.tua, sizeof tcopy);
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t  cc = _popcnt64(tua[i] & BIT_SET_2X);
				vsize[cc].setBit(i);
			}
			nua = 0;
			for (int i = 0; i < 20; i++) if (vsize[i].isNotEmpty()) {
				int j;
				BF128 v = vsize[i];
				while ((j = v.getFirst128()) >= 0) {
					v.clearBit(j);
					register uint64_t U = tcopy[j];
					for (uint32_t k = 0; k < nua; k++) {
						register uint64_t U2 = tua[k];
						if (!(U & (~U2))) { U = 0; break; }
					}
					if (U) tua[nua++] = U;
					if (nua > 40) break;
				}
				if (nua > 40) break;
			}

		}
	}tg2[81], tg3[81], guaw;
	void Init() {
		for (int i = 0; i < 81; i++) {
			tg2[i].Init(0, 0, i); tg3[i].Init(0, 1, i);
		}
	}
	inline void Add2(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg2[i].Add(u); }
	inline void Add3(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg3[i].Add(u); }

	int IsUa4(int i81) {
		GUA& g = tg2[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) == 2) return 1;
		return 0;
	}
	int IsUamin(int i81) {
		GUA& g = tg3[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) < 5) return 1;
		return 0;
	}
	void SortClean3() {
		for (int i = 0; i < 81; i++)
			if (tg3[i].nua > 1) tg3[i].SortClean();
	}
	void SortClean() {
		for (int i = 0; i < 81; i++)
			if (tg2[i].nua > 1) tg2[i].SortClean();
	}
	int CutG2(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg2[i].nua > lim)tg2[i].nua = lim;
			n += tg2[i].nua;
		}
		return n;
	}
	int CutG3(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg3[i].nua > lim)tg3[i].nua = lim;
			n += tg3[i].nua;
		}
		return n;
	}
}guah;

void G17B::UaCollector() {
	zh2gxn.SetupFsol(grid0);
	zh2b[0].InitBands12(grid0);
	zh2gxn.InitKnown(extuas.t2, &extuas.nt2);
	tuasb12.nua = 0;
	FirstUasCollect();// uas 2 digits bands and stacks
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect6();
	Guas2Collect();
	UasCollect4box();// 4 box in bands 1+2
	zh2b[0].InitBands12(grid0);// restore 
	//________________________________ insert stack known uas 
	STD_B3 bw;
	BANDMINLEX::PERM pp;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		int grid_diag[81];
		char ze[81];
		for (int i = 0; i < 81; i++)grid_diag[C_transpose_d[i]] = grid0[i];
		for (int i = 0; i < 81; i++)  ze[i] = grid_diag[i] + '1';
		for (int istack = 0; istack < 3; istack++) {
			int ir = bandminlex.Getmin(&grid_diag[27 * istack], &pp, 0);
			if (ir < 0) {//would be bug} 
				cout << "look for stack bad return" << endl;
				return;
			}
			bw.InitBand2_3(pp.i416, &ze[27 * istack], pp, istack);
			//bw.PrintStatus();
			for (uint32_t i = 0; i < bw.nua; i++) {
				BF128 wd; wd.SetAll_0();// 3x mode to build guam
				uint32_t  zed = 0, ua = bw.tua[i]&BIT_SET_27,ccua=_popcnt32(ua);
				for (uint32_t j = 0; j < 27; j++) {
					if (ua & (1 << j)) {
						int cell = C_transpose_d[j + 27 * istack];
						wd.Set_c(cell);
					}
				}
				int cc3 =  _popcnt32  (wd.bf.u32[2]);
				//cout << Char27out(zed) << endl;
				//cout << Char2Xout(wd.bf.u64[0]) << " ";
				//cout << Char27out(wd.bf.u32[2])<<" "<<wd.Count()<<" "<<cc3 << endl;
				if (cc3 > 3)mb3.Addguam(wd);
				else if(cc3==2  &&  ccua>10) { // must build g2
					int i81=mb3.GetI81_2(wd.bf.u32[2]);
					guah.tg2[i81].AddIf(wd.bf.u64[0]);
				}
				else if (cc3 == 3  &&  ccua>9) { // must build g3
					int i81 = mb3.GetI81_3(wd.bf.u32[2]);
					guah.tg3[i81].AddIf(wd.bf.u64[0]);
				}

			}
		}
		//cout << "band summary ib3 =" << ib3 << "\t";
		//mb3.DumpGuam();


	}
}

inline void G17B::Adduab12( uint32_t digs, uint32_t nd) {
	uint64_t* t = zh2gxn.tua;
	uint32_t n = zh2gxn.nua,n2=0;
	BF128 tw[100];
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t w = t[i], cc = _popcnt64(w);
		if (cc < 25) {	
			BF128 x; x.bf.u64[0] =  w; x.bf.u64[1] = cc;	
			tw[n2++] = x;
		}
	}
	if (n2 > 1) {// sort to check redundancy
		for (uint32_t i = 0; i < n2 - 1; i++) {
			for (uint32_t j=i+1; j < n2 ; j++) {
				if (tw[i].bf.u64[1] > tw[j].bf.u64[1]) {
					BF128 x = tw[i]; tw[i] = tw[j]; tw[j] = x;
				}
			}
		}
	}
	for (uint32_t i = 0; i < n2; i++) {
		BF128 x = tw[i];
		register uint64_t cc = x.bf.u64[1],
			w= x.bf.u64[0],nw=~w;

		if (i) {
			for (uint32_t j = 0; j < i; j++) {
				BF128 y = tw[j];
				register uint64_t ycc = y.bf.u64[1];
				if (ycc == cc) break;
				if(! (y.bf.u64[0] & nw)) { cc = 0; break; }
			}
			if (!cc)continue;// subset found
		}
		if (cc > 15) cc = 15;
		w |= (cc << 59);
		if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w, digs, nd);
	}

}
void G17B::FirstUasCollect() {// produce all uas 2/3 digits
	struct TUAS81 {// used to collect uas 2 digits
		BF128 	tall[200];// first collect
		uint32_t ntall;// , ntold;
		int Add(BF128 w, uint32_t floor) {
			uint32_t cc = w.Count96(), nfloor = ~floor;
			for (uint32_t iua = 0; iua < ntall; iua++) {
				BF128 wt = tall[iua];
				uint32_t cct = wt.bf.u32[3] >> 16;
				uint32_t floort = wt.bf.u32[3] & 0777;
				wt.bf.u32[3] = 0;
				if (cct <= cc)continue; 
				// insert here and check super set
				for (uint32_t jua = ntall; jua > iua; jua--)
					tall[jua] = tall[jua - 1];
				tall[iua] = w;// new inserted
				tall[iua].bf.u32[3] = floor | (cc << 16);
				ntall++;
				return 2;
			}
			w.bf.u32[3] = floor | (cc << 16);
			tall[ntall++] = w;// new added
			return 1;
		}
	}tuas81;
	for (int i = 0; i < 36; i++) {
		uint32_t myf = floors_2d[i];
		zh2_2[0].GoZ2D(myf);
		if (zh2gxn.nua) {
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
			if (cc > 15)cc = 15;
			w |= cc << 59;
			tuasb12.AddInit(w, myf, 2);
			}
		}
	}
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		zhgxn.SetupFsol(grid0);
		for (int i = 0; i < 36; i++) {
			uint32_t myf = floors_2d[i];
			zhou2[0].GoZ2(myf);
			if (zhgxn.nua == 1) {// this is a 18 cells 12+6
				BF128 w = zhgxn.tua[0];
				register uint32_t C = w.bf.u32[2];
				uint32_t  ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777);
				if (ncol != 4)mb3.Addguam(w);
				else {// this is a 6 clues gua2
				}
			}
			else {
				for (uint32_t j = 0; j < zhgxn.nua; j++) {
					BF128 w = zhgxn.tua[j];
					int cc = w.Count();
					int aig = 1;// check for subsets
					if (cc == 18)continue; // has a subset
					for (uint32_t k = 0; k < zhgxn.nua; k++) if (k != j) {
						BF128 w2 = zhgxn.tua[k];
						if (w2.isSubsetOf(w)) { aig = 0; break; }
					}
					if (aig) {
						uint64_t w1 = w.bf.u64[0], cc = _popcnt64(w1);
						register uint32_t C = w.bf.u32[2];
						uint32_t cc2 = _popcnt32(C),
							ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777),
							nplug = 2 * ncol - cc2;
						if ((!cc) || (!cc2)) continue;//ua bands 1+2 or ua band 3
						if (nplug > 2) { mb3.Addguam(w); continue; }
					}
				}
			}		
		}
	}
}
void G17B::SecondUasCollect() {// collect 345 digits in bands 1+2
	for (int i = 0; i < 84; i++) {// look for 3 digits
		uint32_t myf = floors_3d[i];
		int ir = zh2_3[0].GoZ3(myf);// cells unsolved count
		if (ir < 6) continue;// minimum for a fresh ua 3 digits
		uint64_t F = zh2_3[0].cells_unsolved.bf.u64;
		extuas.GetSortA(F);
		zh2_3[0].DoZ3Go();
		if (zh2gxn.nua) Adduab12(myf,3);
	}
	for (int i = 0; i < 126; i++) {// look for 4 digits
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		extuas.GetSortA(F);
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12(floors_4d[i],4);
	}
	for (int i = 0; i < 126; i++) {// look for 5 digits
		uint32_t myf = 0777 ^ floors_4d[i];
		int ir = zh2_5[0].GoZ5(myf);
		uint64_t F = zh2gxn.unsolved_field;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		zh2_5[0].DoZ5Go();
		if (zh2gxn.nua) Adduab12(myf, 5);
	}
}
void G17B::UasCollect6() {
	for (int i = 0; i < 84; i++) {
		uint32_t ass= floors_3d[i], myf = 0777^ass;
		uint64_t bf = zh2gxn.Getsol(ass);
		uint64_t F = bf^BIT_SET_2X;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		//extuas.GetSort(F);
		zh2b[0].InitBf(bf);
		zh2b[0].Dob();
		if (zh2gxn.nua) Adduab12(myf, 6);
	}
}
void G17B::UasCollect4box() {
	//cout << "4 box collect" << endl;
	//___________________ box 1245
	{
		uint64_t stack = units3xBM[5].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		zh2b[0].InitBf(stack);
	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			if (tuasb12.nua < 2550)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	//___________________ box 1346
	{
		uint64_t stack = units3xBM[4].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua,b4);
		zh2b[0].InitBf(stack);
	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	//___________________ box 2356
		{
		uint64_t stack = units3xBM[5].u64[0], b4 = BIT_SET_2X & (~stack);
		//extuas.GetSort(b4);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
}

struct TTT1 {//tuab12;tua <=12 per size, then per floor
#define NB 3
	uint64_t tua0[NB*128], tua[NB * 128];
	uint32_t nua, tdigs0[NB * 128], tdigs[NB * 128],nbl;
	BF128 vs [9][NB] ;// vectors size 4/12
	struct VF {
		BF128 v2[36][NB], v3[84][NB], v4[126][NB], v5[126][NB], v6[84][NB];
	}vf;
	void Init(uint64_t* tu, uint32_t* td, uint32_t n) {
		//cout << "uas <=12 out of n="<<n << endl;
		memset(vs, 0, sizeof vs);
		nua = 0;
		for (uint32_t i = 0; i < n; i++) {
			register uint64_t U = tu[i],cc=U>>59;
			if (cc>=4 && cc <= 12) {
				tua0[nua] = U & BIT_SET_2X;
				tdigs0[nua] = td[i];
				uint32_t bloc = nua / 128, ir = nua - 128 * bloc;
				vs[cc-4][bloc].setBit(ir);
				nua++;
				if (nua >= NB * 128) break;// overflow control
			}
		}
		//cout << "uas <=12 got nua=" << nua << endl;
		// copy per size 
		nbl = (nua + 63) / 64;
		uint32_t x, n2 = 0;
		for (int i = 0; i < 9; i++) {// size 4/12
			uint64_t* v = vs[i][0].bf.u64;
			for (uint32_t j = 0; j < nbl; j++) {// uas for this size
				uint64_t V = v[j];
				while (bitscanforward64(x, V)) {
					V ^= (uint64_t)1 << x;
					register uint32_t ii = x + 64 * j;
					tua[n2] = tua0[ii]  ;
					tdigs[n2++] = tdigs0[ii];
				}
			}
		}
		//DumpUas();
		// build vectors per floor
		memset(&vf, 0, sizeof vf);
		for (int i = 0; i < 36; i++) {// ______________2 digits
			register uint32_t fl = floors_2d[i];
			BF128* pv = vf.v2[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc =j/ 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
		}
		for (int i = 0; i < 84; i++) {// ______________3 digits
			register uint32_t fl = floors_3d[i],fln=~fl;
			BF128* pv = vf.v3[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 36; i++) {// subset___2 digits
				uint32_t fl2 = floors_2d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v2[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________4 digits
			register uint32_t fl = floors_4d[i], fln = ~fl;
			BF128* pv = vf.v4[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 84; i++) {// subset___3 digits
				uint32_t fl2 = floors_3d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v3[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________5 digits
			register uint32_t fl = 0777^floors_4d[i], fln = ~fl;
			BF128* pv = vf.v5[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___4 digits
				uint32_t fl2 = floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v4[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 84; i++) {// ______________6 digits
			register uint32_t fl = 0777 ^ floors_3d[i], fln = ~fl;
			BF128* pv = vf.v6[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___5 digits
				uint32_t fl2 = 0777^floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v5[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		//DumpV();
	}
	int Load(uint64_t* td, int iv, int ix) {
		BF128* v;
		switch (iv) {
		case 2: v = vf.v2[ix]; break;
		case 3: v = vf.v3[ix]; break;
		case 4: v = vf.v4[ix]; break;
		case 5: v = vf.v5[ix]; break;
		case 6: v = vf.v6[ix]; break;
		default: return 0;
		}
		int n = 0,x;
		uint64_t* pv = v->bf.u64;
		for (uint32_t i = 0; i < nbl; i++) {
			uint64_t V = pv[i];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint32_t ii = x + 64 * i;
				td[n++] = tua[ii];
			}
		}
		return n;
	}
	void DumpUas() {
		cout << "uas <=12 sorted" << endl;
		if (nua < 256) for (uint32_t i = 0; i < nua; i++) {
			cout << Char2Xout(tua[i]) << "|";
			cout  << Char9out(tdigs[i]) << " i="<<i << endl;
		}
	}
	void DumpV() {
		cout << " dump v nbl=" << nbl << endl;
		cout << "v2" << endl;
		for (int i = 0; i < 36; i++) {
			uint64_t *pv=vf.v2[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_2d[i]) <<" i="<<i << endl;

		}
		cout << "v3" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v3[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_3d[i]) << " i=" << i << endl;

		}
		cout << "v4" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v4[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v5" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v5[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v6" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v6[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_3d[i]) << " i=" << i << endl;

		}

	}
	//36;84;126;126 for 2;3;4;5; 

}ttt1;

void G17B::Guas2Collect() {
	ttt1.Init(tuasb12.tua, tuasb12.tdigs, tuasb12.nua);
	guah.Init();
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {// initial gua2s 4 cells
		SGUA2 w = tsgua2[i81];
		uint64_t bf = zh2b_g.fd_sols[0][w.dig1].bf.u64;
		bf |= zh2b_g.fd_sols[0][w.dig2].bf.u64;
		uint64_t bf2 = units3xBM[9 + w.col1].u64[0];
		bf2 |= units3xBM[9 + w.col2].u64[0];
		bf &= bf2;
		int xcell1, xcell2;
		bitscanforward64(xcell1, bf);
		bitscanreverse64(xcell2, bf);
		int cell1=From_128_To_81[xcell1], cell2= From_128_To_81[xcell2];
		int r1 = C_row[cell1], r2 = C_row[cell2];
		if (r1 == r2) {
			guah.Add2(bf, i81);
		}
	}
	Guas2CollectG2();
	Guas2CollectG3();

}
void G17B::Guas2CollectG2() {
	//___ find guas2 2 digits
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		int i = 0;
		for (i; i < 36; i++)if (floors_2d[i] == w.digs) break;
		extuas.nt2 = ttt1.Load(extuas.t2, 2, i);
		int ir = zh2_2[0].GoZ2G2(w.digs, w.col1, w.dig1, w.col2, w.dig2);
		if (ir < 0)continue; // locked
		if (ir == 1) {//	solved)
			uint64_t U = zh2gxn.tua[0], cc = _popcnt64(U);
			guah.Add2(U, i81);
		}
		uint64_t F = zh2gxn.unsolved_field;
		zh2_2[0].DoZ2Go();
		if (zh2gxn.nua)
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				uint64_t U = zh2gxn.tua[i];
				guah.Add2(U, i81);
			}
	}
	guah.SortClean();
	//___ find guas2 3 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2=gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			int ir = zh2_3[0].GoZ3G2(fl, w.col1, w.dig1, w.col2, w.dig2),ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 3 digits
	//___ find guas2 4 digits
	guah.SortClean();

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl = floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			int ir = zh2_4[0].GoZ4G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_4[0].DoZ4Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 4 digits

	guah.SortClean();

	//___ find guas2 5 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl =0777^ floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			//if (diag2 == 2)extuas.Dump();
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			int ir = zh2_5[0].GoZ5G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_5[0].DoZ5Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 5 digits
	guah.SortClean();
	ng2=guah.CutG2(30);
}
void G17B::Guas2CollectG3() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs,triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			int i = 0;
			for (i; i < 84; i++)if (floors_3d[i] == w.digs) break;
			extuas.nt2 = ttt1.Load(extuas.t2, 3, i);
			// find UAs 3 digits one 3 digits 
			int ir= zh2_3[0].GoZ3G3(w.digs, rgangcols ),ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add3(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add3(U, i81);
				}
		}
	}
	Guas2CollectG3_4d();
}
void G17B::Guas2CollectG3_4d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do
		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 4, i);
				int ir = zh2_4[0].GoZ4G3(fl, rgangcols), ir2;
				if (ir < 0)continue; // locked
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) guah.Add3(U, i81);						
					}
					continue;
				}
				ir2 = zh2_4[0].DoZ4Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) guah.Add3(U, i81);						
					}
			}
		}
	}
	guah.SortClean3();
	Guas2CollectG3_5d();
}
void G17B::Guas2CollectG3_5d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do
		if (gt.nua >= 10) continue;

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = 0777^floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 5, i);
				int ir = zh2_5[0].GoZ5G3(fl, rgangcols), ir2;
				if (ir < 0)continue; // locked
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) {
							guah.Add3(U, i81);
						}
					}
					continue;
				}
				ir2 = zh2_5[0].DoZ5Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) {
							guah.Add3(U, i81);
						}
					}
			}
		}
	}
	guah.SortClean3();
	ng3 = guah.CutG3(20);
}

//__________ end of initial uas harvest switch to 54 mode

void G17B::StartAfterUasHarvest() {
	guah54.Build();
	t54b12.Build_ta128(tuasb12.tua, tuasb12.nua);
	Expand_03();
}

void T54B12::Build_ta128(uint64_t* t, uint32_t n) {
	// here insert all missing uas one band (1 or 2)
	for (int i = 0; i < 2; i++) {
		STD_B416& wb = (i) ? myband2 : myband1;
		BF64 wt; wt.bf.u64 = 0;
		uint32_t* tu = wb.tua, nu = wb.nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			wt.bf.u32[i] = u;
			if (cc > 12)if(n< UA12SIZE)t[n++] = wt.bf.u64;
		}
	}	
	InitA();
	//BF128 vsize[25][30];
	uint32_t nbl64 = (n + 63) >> 6, x;
	memset(vsize, 0, sizeof vsize);
	for (uint32_t i = 0; i < n; i++) {
		int bloc = i >> 7, ir = i - (bloc << 7);
		uint64_t cc = t[i] >> 59;
		if (cc > 24) continue;//safety
		vsize[cc][bloc].setBit(ir);
	}
	//uint64_t tw[UA12SIZE]; // to check redundancy moved to t54b12
	uint32_t ntw = 0;
	for (int i1 = 4; i1 < 25; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {// switch to 54 mode here
				V ^= (uint64_t)1 << x;
				register uint64_t R = t[x + (i2 << 6)];
				R = (R & BIT_SET_27) | ((R & BIT_SET_B2) >> 5);
				if (1) {// check redundancy and subsets
					for (uint32_t i = 0; i < ntw; i++) {
						if (!(R & (~tw[i]))) {
							//cout << Char54out(R) << " erased" << endl;;
							//cout << Char54out(tw[i]) <<" i="<<i << endl;
							R = 0; break;						}
					}
					if (R)tw[ntw++] = R;
					else continue;
				}
				AddA(R);
			}
		}
	}
}
int T54B12::Build_tb128() {
	int tc[3], ntc = 0;
	{// // build table of clues 
		int cell;
		register uint64_t U = spb_0_15[3].all_previous_cells;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}
	//BF128 tvw[20];
	uint32_t lastbloc = t54b12.nablocs;
	tvw[0] = spb_0_15[3].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = ta128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	//BF128 vsize[19][UABNBLOCS];
	memset(vsize, 0, sizeof vsize);
	//uint64_t tw[128 * UABNBLOCS];
	uint32_t ntw = 0;
	{
		register uint64_t Ac = spb_0_15[3].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = ta128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UABNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UABNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors
	InitB();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				AddB(tw[x + (i2 << 6)]);
			}
		}
	}
	if (nb128 < 128) nb128 = 128; //forcing add outside first bloc
	return 0;
}
int T54B12::Build_tc128() {
	int tc[3], ntc = 0;
	{// // build table of clues 
		int cell;
		register uint64_t U = spb_0_15[7].all_previous_cells
			& (~spb_0_15[3].all_previous_cells);// fresh clues 
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}
	//BF128 tvw[UABNBLOCS];
	uint32_t lastbloc = t54b12.nbblocs;
	tvw[0] = spb_0_15[7].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = tb128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	//BF128 vsize[19][UACNBLOCS];
	memset(vsize, 0, sizeof vsize);
	//uint64_t tw[128 * UACNBLOCS];
	uint32_t ntw = 0;
	{
		register uint64_t Ac = spb_0_15[7].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tb128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UACNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UACNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitC();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in tc128[0]
				if (t54b12.IsNotRedundant(U))
					AddC(U);
			}
		}
	}
	if (nc128 < 128) nc128 = 128; //forcing add outside first bloc
	return 0;
}
int T54B12::Build_td128() {
	int tc[3], ntc = 0;
	{// // build table of clues  ater 3+3+3
		int cell;// 0 123 4 567  8  9 10 11    12 13 14 15  
		register uint64_t U = g17b.myb12_9
			& (~spb_0_15[7].all_previous_cells);// fresh clues 
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			tc[ntc++] = cell;
		}
	}
	//BF128 tvw[UACNBLOCS];
	uint32_t lastbloc = t54b12.ncblocs;
	for (uint32_t i = 0; i <= lastbloc; i++) {
		TUVECT& vv = tc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[tc[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	//BF128 vsize[19][UADNBLOCS];
	memset(vsize, 0, sizeof vsize);
	//uint64_t tw[128 * UADNBLOCS];
	uint32_t ntw = 0;
	{
		register uint64_t Ac = g17b.myac_9;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tc128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UADNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UADNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitD();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in tc128[0]
				if (t54b12.IsNotRedundantD(U))
					AddD(U);
			}
		}
	}
	return 0;
}

void GUAH54::Build() {// cut to 30 switch to 54 killer

	uint64_t* pbuf = gbuf;
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			GUA54& gd = tg2[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 50;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUAH::GUA& g0 = guah.tg3[i81];
			GUA54& gd = tg3[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
}
void GUAH54::Build2(uint64_t filter, uint64_t active) {
	register uint64_t F = filter,
		A = active;

	uint64_t* pbuf = gbuf;
	uint64_t tw[60];// valid
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUA54& g0 = guah54.tg2[i81];
			GUA54& gd = tg2[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20)cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 20;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUA54& g0 = guah54.tg3[i81];
			GUA54& gd = tg3[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20) cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
}
void GUAH54::Build9(uint64_t filter, uint64_t active) {
	register uint64_t F = filter,		A = active;
	uint64_t* pbuf = gbuf;
	uint64_t tw[60];// valid
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUA54& g0 = guah54_2.tg2[i81];
			GUA54& gd = tg2[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				//pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20)cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 20;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUA54& g0 = guah54_2.tg3[i81];
			GUA54& gd = tg3[i81];
			if (g0.killer & F) {
				pbuf += gd.nuamax;// lock the storing place
				continue; // killed or not valid
			}
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf += gd.nuamax;// lock the storing place
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				if (cc > 20) cc = 20;
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g3 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i <= 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
				if (gd.nua >= 30) break;
			}
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
}


BF128 GUAH54::GetG2(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g2;
	g2.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg2[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) {	g2.setBit(i81);	continue;}
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g2.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g2;
}
BF128 GUAH54::GetG3(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g3;
	g3.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg3[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) { g3.setBit(i81);	continue; }
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g3.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g3;
}

///________ start expand uas bands 12
/// target 10 clues 3+3+4
/// target 11 clues 3+3+5
/// target 12 clues 3+3+3+3
/// last step storing potential valid bands 1+2 of the step
/// end of early steps skrinking and restructuring uas bands 1+2
/// and building a reduced GUAs table 

void G17B::Expand_03() {
	SPB03* s, * sn;
	if (aigstop) return;
	T54B12::TUVECT& tuv128 = t54b12.ta128[0];
	uint64_t* twu = tuv128.t;
	s = spb_0_15;		memset(s, 0, sizeof spb_0_15[0]);
	s->active_cells = maskLSB[54].u64[0];
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done

	//_______ start search 3 first clues
next:
	// catch and apply cell in bitfields
	int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= spb_0_15)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tuv128.vc[cell];
	if (sn->ncl == 3) {// 3 cellsfirst step
		p_cpt2g[3]++;
		if (t54b12.Build_tb128()) goto next;
		Expand_46();		
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	if (ir < 0) return;//never
	uint64_t Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch unlikely
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}
void G17B::Expand_46() {
	if (aigstop) return;
	SPB03*sl= &spb_0_15[4] ,* s=sl, * sn;
	T54B12::TUVECT& tu128 = t54b12.tb128[0];
	uint64_t* twu = tu128.t;
	*s = spb_0_15[3];	// duplicate 3 for new vector
	s->possible_cells = twu[0];
	s->v = tu128.v0;// initial nothing done
	//_______ start search clues 4-6
next:	// catch and apply cell in bitfields
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= sl)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->all_previous_cells |= bit;
	sn->v &= tu128.vc[cell];
	if (sn->ncl == 6) {// 6 cells
		p_cpt2g[4]++;
		if (t54b12.Build_tc128()) goto next;
		GoExpand_7_10();
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	uint64_t Ru;
	if (ir < 0) {//never valid not a critical path
		uint32_t tc[8], ntc = 0;
		{// // build a safe table of clues 
			int cell;
			register uint64_t U = sn->all_previous_cells;// fresh clues 
			while (bitscanforward64(cell, U)) {
				U ^= (uint64_t)1 << cell;
				tc[ntc++] = cell;
			}
		}
		if (t54b12.nbblocs) {// more uas to check
			for (uint32_t i = 1; i <= t54b12.nbblocs; i++) {
				T54B12::TUVECT& vv = t54b12.tb128[i];
				BF128 v = vv.v0, * vc = vv.vc;
				for (uint32_t ic = 0; ic < ntc; ic++)
					v &= vc[tc[ic]];
				if (v.isNotEmpty()) {
					int ir2 = v.getFirst128();
					uint64_t Ru = vv.t[ir2] & ac;
					if (!Ru)goto next;//dead branch
					sn->possible_cells = Ru;
					s++; // switch to next spot
					goto next;
				}
			}

		}
		if (zh2b[1].IsValid(sn->all_previous_cells)) {
			uint32_t i = zh2gxn.nua - 1;
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			t54b12.AddA(ua54);
			t54b12.AddB(ua54);
			Ru = ua;
		}
		else {
			cout << "bug exp 4-7 lower 7" << endl;
			aigstop = 1; return;
		}

	}
	else  Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}

int G17B::IsValid7pbf(SPB03* sn) {
	if (zh2b[1].IsValid(sn->all_previous_cells)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 21) {
				t54b12.AddA(ua54);
				t54b12.AddB(ua54);
				t54b12.AddC(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
inline int G17B::GetNextCell(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->all_previous_cells |= bit;
	sn->v &= t54b12.tc128[0].vc[cell];
	return 0;
}
inline void G17B::GetNextUa(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v.bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v.bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir + 64];
	}
}
inline void G17B::GetNextUaAdd(SPB03* sn) {
	if (t54b12.nc128 <= 128) return;
	// more uas to check
	for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
		T54B12::TUVECT& vv = t54b12.tc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 6; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return;
		}
	}
}
inline int G17B::GetLastAndUa(SPB03* sn, int d) {
	int aig = 0;
	register uint64_t  V, And = ~0;
	register uint32_t ir;
	if ((V = sn->v.bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir];
		}
	}
	if (!And) return 1;
	if ((V = sn->v.bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (t54b12.ncblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			T54B12::TUVECT& vv = t54b12.tc128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 6; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}

/// last step in bands 1+2 expansion critical code
/// expand builds tables of potential solutions 
/// one with target number of clues stored 
///		 bit field ;and of all items; number of clues/bands stacks)
/// one per number of clues below the target 
///      storing also the active cells downstream

void G17B::GoExpand_7_10() {// here only 18 clues search 
	if (op.t18) {
		if ( op.p1) { GoExpand_7_11(); return; }
		else { GoExpand_7_12(); return; }
	}
	return; 
}

/// 7_11 is for 18 pass1 
/// stack limit bands limit is 6 in pass 2

void G17B::Go_10_11_18() {// 10 clues limit 11 clues 
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		BF128 ww = tbelow10[iv];
		cb3.ncl = 10;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 10);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		// can be here 288 or 828
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) 
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		{
			register uint64_t nb1 = cb3.cbs.b[0], nb2 = cb3.cbs.b[1];
			if (nb1 > 7 || nb2 > 7) continue; 
			if (nb1 == 7) Ac &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) Ac &= BIT_SET_27;
		}
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 11;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);

		}
	}
}
void G17B::Go_9_11_18() {// 9 clues limit 11 clues 
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		cb3.ncl = 9;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 9);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 7+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 10;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			// could be 288 828
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 7 || nb2 > 7) continue;
				if (nb1 == 7) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 7) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 11;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}
		}
	}

}
void G17B::Go_8_11_18() {// 8 clues limit 11 clues 
	//cout << " entry 8 clues for 11 clues" << endl;
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		cb3.ncl = 8;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 8);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 6+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 9;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 10;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				// could be 288 828
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 7 || nb2 > 7) continue;
					if (nb1 == 7) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 7) Ac3 &= BIT_SET_27;
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 11;
					cb3n3.cbs.Add(cell3);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
				}
			}
		}
	}

}
void G17B::Go_7_11_18() {// 7 clues limit 11 clues 
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[0]; iv++) {
		BF128 ww = tbelow7[iv];
		cb3.ncl = 7;
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 7);
		// try direct
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		p_cpt2g[7]++;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
			genb12.bands3[ib3].Go(cb3);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;// maxi is 5+2 possible 8+2
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.ncl = 8;
			cb3n.cbs.Add(cell);
			cb3n.g2t = guah54_2.GetG2(myb12);
			cb3n.g3t = guah54_2.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.ncl = 9;
				cb3n2.cbs.Add(cell2);
				cb3n2.g2t = guah54_2.GetG2(myb12);
				cb3n2.g3t = guah54_2.GetG3(myb12);
				// could be 288 828
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 10;
					cb3n3.cbs.Add(cell3);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
					// try now a fourth  clue in bands 1+2
					uint64_t Ac4 = Ac3;
					{
						register uint64_t nb1 = cb3n3.cbs.b[0], nb2 = cb3n3.cbs.b[1];
						if (nb1 > 7 || nb2 > 7) continue;
						if (nb1 == 7) Ac4 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 7) Ac4 &= BIT_SET_27;
					}
					int cell4;
					while (bitscanforward64(cell4, Ac4)) {
						CALLBAND3 cb3n4 = cb3n3;
						uint64_t bit4 = (uint64_t)1 << cell4;
						Ac4 ^= bit4; //clear bit
						cb3n4.bf12 |= bit4;
						myb12 = cb3n4.bf12;
						cb3n4.ncl = 11;
						cb3n4.cbs.Add(cell4);
						cb3n4.g2t = guah54_2.GetG2(myb12);
						cb3n4.g3t = guah54_2.GetG3(myb12);
						for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
							genb12.bands3[ib3].Go(cb3n4);
					}
				}
			}
		}
	}	
}

int G17B::Expand_7_11() {
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	uint32_t nfull = 0, nbelow = 0;
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 7 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	//_______ start search clues 7_11
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	
	else {	return (nfull | nbelow);	}
	sn = s + 1;
	ua_ret7p = 0;
	
	if (sn->ncl == 11) {// 11 cells
		tfull[ntbelow[5]++] = sn->all_previous_cells;
		tandbelow[5] &= sn->all_previous_cells;
		nfull++;
		goto next;
	}
	
	if (sn->ncl == 10) {// last step
		p_cpt2g[5]++;
		if (!GetLastAndUa(sn)) {// can be a valid 10 or unknown ua(s)
			if (IsValid7pbf(sn)) // got uas to use
				ua_ret7p = anduab12;
			else {// valid 10
				p_cpt2g[57]++;
				BF128 w;
				w.bf.u64[0] = sn->all_previous_cells;
				w.bf.u64[1] = sn->active_cells;
				tbelow10[ntbelow[3]++] = w;
				tandbelow[3] &= sn->all_previous_cells;			
				nbelow++;
				goto next;
			}

		}
		{// last not empty direct or after check valid
			if (!ua_ret7p) goto next;
			register uint64_t nb1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
					nb2 = 10 - nb1, P = ua_ret7p & s->active_cells;
			if (nb1 > 7 || nb2 > 7) goto next;
			if (nb1 == 7) P &= ~(uint64_t)BIT_SET_27;
			if (nb2 == 7) P &= BIT_SET_27;
			if (P) {	sn->possible_cells = P; s++;	}
			goto next;
		}

	}
	else {
		if (sn->v.isNotEmpty())GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {// add ua or store below
		if(!IsValid7pbf(sn))  {// valid below
			clean_valid_done = 1;
			int isize = sn->ncl - 7;//table index 0 7 clues
			p_cpt2g[54 + isize]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			switch (isize) {
			case 0: tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= sn->all_previous_cells; break;
			case 1: tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= sn->all_previous_cells; break;
			case 2: tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells; break;
			}
			nbelow++;
			goto next;
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}
void G17B::GoExpand_7_11() {

	SPB03* sn = &spb_0_15[7];
	//nclues = 6;// see why needed
	t54b12.ntm = 0;
	if (!Expand_7_11()) return;
	p_cpt2g[6]++;
	int nbelow = ntbelow[4] + ntbelow[3] + ntbelow[2] + ntbelow[1] + ntbelow[0];
	p_cpt2g[51] += ntbelow[5];
	p_cpt2g[52] += nbelow;
	if (p_cpt2g[61] < ntbelow[5])p_cpt2g[61] = ntbelow[5];
	if (p_cpt2g[62] < ntbelow[3])p_cpt2g[62] = ntbelow[3];

	myandall = tandbelow[0]& tandbelow[1]& tandbelow[2]& tandbelow[3]& tandbelow[4]& tandbelow[5];
	guah54_2.Build2(myandall, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) 
		genb12.bands3[ib3].BuildGuam2(myandall);
	if (ntbelow[0])  Go_7_11_18(); // do 7 clues then more
	if (ntbelow[1])  Go_8_11_18(); // do 8 clues then more
	if (ntbelow[2])  Go_9_11_18(); // do 9 clues then more
	if (ntbelow[3])  Go_10_11_18(); // do 10 clues then more

	if (ntbelow[5] == 1) {// get active g2 g3 from guah54_2 direct
		p_cpt2g[7]++;
		cb3.ncl = 11;
		myb12 = cb3.bf12 = tfull[0];
		cb3.cbs.Init(myb12, 11);
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			genb12.bands3[ib3].Go(cb3);
			if (clean_valid_done == 2)break;
		}
	}
	else {
		for (uint32_t iv = 0; iv < ntbelow[5]; iv++) {
			cb3.ncl = 11;
			myb12 = cb3.bf12 = tfull[iv];
			if (t54b12.NotValid(myb12)) continue;// uas b12 added
			cb3.cbs.Init(myb12, 11);
			cb3.g2t = guah54_2.GetG2(cb3.bf12);
			cb3.g3t = guah54_2.GetG3(cb3.bf12);
			p_cpt2g[7]++;
			clean_valid_done = 0;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				genb12.bands3[ib3].Go(cb3);
				if (clean_valid_done == 2)break;
			}
		}
	}
}

//_________________________ processing 18 clues pass 2 666 666

inline int G17B::GetNextCellD(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= t54b12.td128[0].vc[cell];
	return 0;
}
inline void G17B::GetNextUaD(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v.bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.td128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v.bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.td128[0].t[ir + 64];
	}
}
inline void G17B::GetNextUaAddD(SPB03* sn) {
	if (t54b12.nd128 <= 128) return;
	// more uas to check
	for (uint32_t i = 1; i <= t54b12.ndblocs; i++) {
		T54B12::TUVECT& vv = t54b12.td128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 9; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return;
		}
	}
}
inline int G17B::GetLastAndUaD(SPB03* sn, int d) {
	int aig = 0;
	register uint64_t  V, And = sn->active_cells;
	register uint32_t ir;
	if ((V = sn->v.bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.td128[0].t[ir];
		}
	}
	if (!And) return 1;
	if ((V = sn->v.bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.td128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (t54b12.ndblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= t54b12.ndblocs; i++) {
			T54B12::TUVECT& vv = t54b12.td128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 9; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}

void G17B::Go_11_12() {
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[4]; iv++) {
		BF128 ww = tbelow11[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 11);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		Ac &= (cb3.cbs.NextActive() & cb3.cbs.NextActiveStack());
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			cb3n.ncl = 12;
			cb3n.g2t = guah54_9.GetG2(myb12);
			cb3n.g3t = guah54_9.GetG3(myb12);
			p_cpt2g[7]++;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
				genb12.bands3[ib3].Go(cb3n);
		}
	}
}
void G17B::Go_10_12() {
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[3]; iv++) {
		BF128 ww = tbelow10[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 10);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			Ac2 &= (cb3n.cbs.NextActive() & cb3n.cbs.NextActiveStack());	
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				cb3n2.ncl = 12;
				cb3n2.g2t = guah54_9.GetG2(myb12);
				cb3n2.g3t = guah54_9.GetG3(myb12);
				p_cpt2g[7]++;
				for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
					genb12.bands3[ib3].Go(cb3n2);
			}
		}
	}
}

void G17B::Go_9_12() {
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[2]; iv++) {
		BF128 ww = tbelow9[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 9);
		// try now one more clue in bands 1+2
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					{
						register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
						if (nb1 > 6 || nb2 > 6) continue;
						if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 6) Ac3 &= BIT_SET_27;
						register uint64_t S = stack1_54;
						if (cb3n2.cbs.s[0] == 6) Ac3 &= ~S;
						if (cb3n2.cbs.s[1] == 6) Ac3 &= ~(S << 3);
						if (cb3n2.cbs.s[2] == 6) Ac3 &= ~(S << 6);
					}
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.ncl = 12;
					cb3n3.cbs.Add(cell2);
					cb3n3.g2t = guah54_2.GetG2(myb12);
					cb3n3.g3t = guah54_2.GetG3(myb12);
					p_cpt2g[7]++;
					for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
						genb12.bands3[ib3].Go(cb3n3);
				}
			}
		}
	}
}
void G17B::Go_8_12() {
	t54b12.ntm = 0;
	clean_valid_done = 1;
	for (uint32_t iv = 0; iv < ntbelow[1]; iv++) {
		BF128 ww = tbelow8[iv];
		myb12 = cb3.bf12 = ww.bf.u64[0];
		cb3.cbs.Init(myb12, 8);
		uint64_t Ac = ww.bf.u64[1];
		int cell;
		while (bitscanforward64(cell, Ac)) {
			CALLBAND3 cb3n = cb3;
			uint64_t bit = (uint64_t)1 << cell;
			Ac ^= bit; //clear bit
			cb3n.bf12 |= bit;
			myb12 = cb3n.bf12;
			cb3n.cbs.Add(cell);
			// try now a second clue in bands 1+2
			uint64_t Ac2 = Ac;// others are not active now
			{
				register uint64_t nb1 = cb3n.cbs.b[0], nb2 = cb3n.cbs.b[1];
				if (nb1 > 6 || nb2 > 6) continue;
				if (nb1 == 6) Ac2 &= ~(uint64_t)BIT_SET_27;
				if (nb2 == 6) Ac2 &= BIT_SET_27;
			}
			int cell2;
			while (bitscanforward64(cell2, Ac2)) {
				CALLBAND3 cb3n2 = cb3n;
				uint64_t bit2 = (uint64_t)1 << cell2;
				Ac2 ^= bit2; //clear bit
				cb3n2.bf12 |= bit2;
				myb12 = cb3n2.bf12;
				cb3n2.cbs.Add(cell2);
				// try now a third clue in bands 1+2
				uint64_t Ac3 = Ac2;// others are not active now
				{
					register uint64_t nb1 = cb3n2.cbs.b[0], nb2 = cb3n2.cbs.b[1];
					if (nb1 > 6 || nb2 > 6) continue;
					if (nb1 == 6) Ac3 &= ~(uint64_t)BIT_SET_27;
					if (nb2 == 6) Ac3 &= BIT_SET_27;
				}
				int cell3;
				while (bitscanforward64(cell3, Ac3)) {
					CALLBAND3 cb3n3 = cb3n2;
					uint64_t bit3 = (uint64_t)1 << cell3;
					Ac3 ^= bit3; //clear bit
					cb3n3.bf12 |= bit3;
					myb12 = cb3n3.bf12;
					cb3n3.cbs.Add(cell3);
					// try now a fourth  clue in bands 1+2
					uint64_t Ac4 = Ac3;
					{
						register uint64_t nb1 = cb3n3.cbs.b[0], nb2 = cb3n3.cbs.b[1];
						if (nb1 > 6 || nb2 > 6) continue;
						if (cb3n3.cbs.s[0] > 6 || cb3n3.cbs.s[1] > 6 || cb3n3.cbs.s[2] > 6) continue;
						if (nb1 == 6) Ac4 &= ~(uint64_t)BIT_SET_27;
						if (nb2 == 6) Ac4 &= BIT_SET_27;
						register uint64_t S = stack1_54;
						if (cb3n3.cbs.s[0]== 6)Ac4 &= ~S;
						if( cb3n3.cbs.s[1] == 6 )Ac4 &= ~(S<<3);
						if( cb3n3.cbs.s[2] == 6) Ac4 &= ~(S << 6);
					}
					int cell4;
					while (bitscanforward64(cell4, Ac4)) {
						CALLBAND3 cb3n4 = cb3n3;
						uint64_t bit4 = (uint64_t)1 << cell4;
						Ac4 ^= bit4; //clear bit
						cb3n4.bf12 |= bit4;
						myb12 = cb3n4.bf12;
						cb3n4.ncl = 12;
						cb3n4.cbs.Add(cell4);
						cb3n4.g2t = guah54_2.GetG2(myb12);
						cb3n4.g3t = guah54_2.GetG3(myb12);
						p_cpt2g[7]++;
						for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
							genb12.bands3[ib3].Go(cb3n4);
					}
				}
			}
		}
	}
}



void G17B::Expand_7_9() {
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	tand_7_9 = ~0;
	nt_7_9 = 0;
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done


	//_______ start search clues 7_9
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	else 	return ;
	sn = s + 1;
	{// adjust active to the band constraint
		register uint64_t nb1 = _popcnt64(sn->all_previous_cells & BIT_SET_27),
			nb2 = sn->ncl - nb1, P =  sn->active_cells;
		if (nb1 > 6 || nb2 > 6) goto next;
		if (nb1 == 6) P &= ~(uint64_t)BIT_SET_27;
		if (nb2 == 6) P &= BIT_SET_27;
		sn->active_cells=P;
	}
	ua_ret7p = 0;
	if (sn->v.isNotEmpty())GetNextUa(sn);// first 128
	else GetNextUaAdd(sn);
	if (!ua_ret7p) {// add ua or go below
		if (!IsValid7pbf(sn)) {// valid below
			int isize = sn->ncl - 7;//table index 0 7 clues
			p_cpt2g[54 + isize]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			switch (isize) {
			case 0: tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= sn->all_previous_cells; break;
			case 1: tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= sn->all_previous_cells; break;
			case 2: tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells; break;
			}
			goto next;
		}
	}
	// now a new ua is granted 
	if (sn->ncl == 9) {// 9 cells store it
		BF128 w;
		w.bf.u64[0] = sn->all_previous_cells;
		w.bf.u64[1] = sn->active_cells;
		t_7_9[nt_7_9++] = w;
		tand_7_9 &= sn->all_previous_cells;
		goto next;
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}
void G17B::Expand_10_12() {
	SPB03* sl = &spb_0_15[12], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.td128[0];
	uint64_t* twu = tuv128.t;
	//*s initial done by the caller
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	//_______ start search clues 7_11
next:	// catch and apply cell in bitfields
	if (GetNextCellD(s))	if (--s >= sl)goto next;else  return ;
	sn = s + 1;
	//cout << Char54out(sn->all_previous_cells);
	//sn->cbs.Status(); cout << " ncl=" << sn->ncl<<endl;
	ua_ret7p = 0;
	// adjust active to the band constraint
	if (!(sn->active_cells &= sn->cbs.NextActive())) goto next;
	if (sn->ncl == 11) {// last step
		p_cpt2g[5]++;
		if (!GetLastAndUaD(sn)) {// can be a valid 10 or unknown ua(s)
			if (IsValid7pbf(sn)) {// got uas to use add it to D	
				t54b12.AddD(ua_ret7p);// re use smallest
				ua_ret7p = anduab12;
			}
			else {// valid 11
				p_cpt2g[58]++;
				BF128 w;
				w.bf.u64[0] = sn->all_previous_cells;
				w.bf.u64[1] = sn->active_cells;
				tbelow11[ntbelow[4]++] = w;
				tandbelow[4] &= sn->all_previous_cells;
				goto next;
			}
		}
		{// last not empty direct or after check valid
			register uint64_t  P = ua_ret7p & sn->active_cells;
			P &= sn->cbs.NextActiveStack();
			// expand directly P to 12, these are potential 12 666 666
			int cell;
			register uint64_t pc = sn->all_previous_cells;
			while (bitscanforward64(cell, P)) {
				register uint64_t bit  = (uint64_t)1 << cell;
				P ^= bit;
				bit |= pc;
				tfull[ntbelow[5]++] = bit;
				tandbelow[5] &= bit;
			}
			goto next;
		}
	}
	else {
		if (sn->v.isNotEmpty())GetNextUaD(sn);// first 128
		else GetNextUaAddD(sn);
	}

	if (!ua_ret7p) {// can only be 10 clues
		if (!IsValid7pbf(sn)) {// valid below
			p_cpt2g[56]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			tbelow10[ntbelow[3]++] = w;
			tandbelow[3] &= sn->all_previous_cells;
			goto next;
		}
		else 	t54b12.AddD(ua_ret7p);// re use smallest		
	}
	sn->possible_cells = ua_ret7p & sn->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}

void G17B::GoExpand_7_12() {
	SPB03* sn = &spb_0_15[7];
	t54b12.ntm = 0;
	p_cpt2g[70]++;

	Expand_7_9();

	if (ntbelow[0]) {
		cout << "unexpected 7 clues pass2 stop [4]" << p_cpt2g[4] << endl;
		aigstop = 1;		return;
	}
	p_cpt2g[71]++;
	if (p_cpt2g[72] < nt_7_9)p_cpt2g[72] = nt_7_9;
	myandall = tandbelow[1] & tandbelow[2]& tand_7_9;
	guah54_2.Build2(myandall, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam2(myandall);

	if (ntbelow[1]) Go_8_12();
	if (ntbelow[2]) Go_9_12();
	for (uint32_t i1 = 0; i1 < nt_7_9; i1++)  
		GoExpand_10_12(t_7_9[i1]);
}
void G17B::GoExpand_10_12(BF128 ww){
	myb12_9 = ww.bf.u64[0];
	myac_9 = ww.bf.u64[1];
	CBS cbs;
	if (t54b12.Build_td128()) return;
	cbs.Init(myb12_9, 9);
	spb_0_15[12].Init9(ww, cbs);
	// clear processed 7,8,9 init the count for next expansion step
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	if (t54b12.nd128 < 128)t54b12.nd128 = 128;// force adds outside
	Expand_10_12();
	p_cpt2g[73] += ntbelow[5];
	if (p_cpt2g[74] < ntbelow[5]) p_cpt2g[74] = ntbelow[5];

	if ((!ntbelow[5]) && (!ntbelow[4]) && (!ntbelow[3])) return;

	myandall_9 = tandbelow[3] & tandbelow[4] & tandbelow[5];
	guah54_9.Build9(myandall_9, myac_9);

	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].BuildGuam9(myandall);

	if (ntbelow[4])Go_11_12();
	if (ntbelow[3])Go_10_12();


	if (ntbelow[5] == 1) {// get active g2 g3 from guah54_2 direct
		p_cpt2g[7]++;
		cb3.ncl = 12;
		myb12 = cb3.bf12 = tfull[0];
		cb3.cbs.Init(myb12, 12);
		cb3.g2t = guah54_9.GetG2(cb3.bf12);
		cb3.g3t = guah54_9.GetG3(cb3.bf12);
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			genb12.bands3[ib3].Go(cb3);
			if (clean_valid_done == 2)break;
		}

	}
	else {
		for (uint32_t iv = 0; iv < ntbelow[5]; iv++) {
			cb3.ncl = 12;
			myb12 = cb3.bf12 = tfull[iv];
			if (t54b12.NotValid(myb12)) continue;// uas b12 added
			cb3.cbs.Init(myb12, 12);
			cb3.g2t = guah54_9.GetG2(cb3.bf12);
			cb3.g3t = guah54_9.GetG3(cb3.bf12);
			p_cpt2g[7]++;
			clean_valid_done = 0;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				genb12.bands3[ib3].Go(cb3);
				if (clean_valid_done == 2)break;
			}
		}
	}
}


//____________processing band 3 to the end


int tstack27[3] = { 07007007 ,070070070,0700700700 };
struct CRITB3 {
	uint32_t minix[4],// triplet bf1 bf2 bf3  
		critbf, pairs27, mincount,
		critstack,
		assigned, active,
		ncl, nb3, nmiss;
	BF64 stackf, stackb12, stackmin;
	inline void Init(int ncb12,  CBS& cbsx) {
		memset(this, 0, sizeof(*this));
		ncl = ncb12;
		nb3 = 18 - ncl;
		memcpy(&stackb12, cbsx.s, sizeof stackb12);
	}
	inline int Addone(uint32_t i27) {// back 1 if not possible
		int bit27 = 1 << i27, stack = C_stack[i27];
		int imini = i27 / 3, bitmini = 1 << imini, bitstack = 1 << stack;
		assigned |= bit27;
		if (!(bit27 & critbf)) {// clue added outfield
			if (!nmiss)return 1;// not possible
			if (critstack & bitstack)return 1;// not possible  
			nmiss--;
			if (!nmiss) {
				critstack = 7;
				active &= critbf;// no more outfield
			}
			else {
				stackf.bf.u16[stack]++;
				if (stackf.bf.u16[stack] >= nb3) {
					critstack |= bitstack;
					// limit the stack to bit field
					register int mask =~( 07007007 << (3 * stack));
					mask |= critbf;
					active  &= mask;
				}
			}
			return 0;;
		}
		// now add in field within mincount
		if (minix[3] & bitmini) {// 2 clues expected
			critbf ^= bit27; minix[3] ^= bitmini;
			if (!nmiss)active &= critbf;
			return 0;
		}
		if (minix[2] & bitmini){
			if (pairs27 & bit27) {// 2 clues not common clue one more clue
				if (critstack & bitstack)return 1;// not possible
				critbf ^= bit27;	minix[2] ^= bitmini;
				nmiss--;
				if (!nmiss) {
					critstack = 7;
					active &= critbf;// no more outfield
				}
				else {
					stackf.bf.u16[stack]++;
					if (stackf.bf.u16[stack] == nb3) {
						critstack |= bitstack;
						// limit the stack to bit field
						register int mask = ~(07007007 << (3 * stack));
						mask |= critbf;
						active &= mask;
					}
				}
				return 0;
			}
		}
		register int mask = ~(7 << (3 * imini));// clear minirow
		critbf &= mask;		
		if( !nmiss )active &= critbf;
		else if(critstack & bitstack)active &= mask;
		if (minix[1] & bitmini) minix[1] ^= bitmini; // one clue expected 		
		else if (minix[2] & bitmini)minix[2] ^= bitmini;//common cell in 2 pairs 
		else if (minix[0] & bitmini)minix[0] ^= bitmini;// triplet
		return 0;
	}
	inline int AddAssign(uint32_t bfa) {// back 1 if not possible
		if (assigned & bfa)return 1; //should never be
		active &= ~bfa; // minimum is to kill new assign
		register int i27, X = bfa;
		while (bitscanforward(i27, X)) {
			X ^= 1 << i27;// clear bit
			if (Addone(i27))return 1;
		}
		return 0;
	}
	void AssignBf2(int bf) {// all or part of the bf2
		minix[2] &= ~bf;// clean the bf2 
		int imini;
		while (bitscanforward(imini, bf)) {
			int  mask = 7 << (3 * imini), m27 = pairs27 & mask;
			bf ^= 1 << imini;// clean the bf
			critbf &= ~mask;// clean the mini row as crit field
			active &= ~mask;// no more clue in this minirow
			pairs27 ^= m27;// clean the 27 pairs to assign
			assigned |= mask ^ m27;// assign the third cell of the mini row
		}
	}
	inline void AssignCritical() {
		active = critbf;
		if (critstack == 7) {
			if (minix[2])  AssignBf2(minix[2]);
			return;
		}
		else if (critstack) {
			for (int i = 0, bit = 1, mask = 07007007; i < 3; i++, bit <<= 1, mask <<= 3)
				if (!(critstack & bit)) active |= mask;
		}
		else active = BIT_SET_27;
		if ((!critstack) || (!minix[2])) return;
		for (int i = 0, bit = 1, mask = 0111; i < 3; i++, bit <<= 1, mask <<= 1) {
			if (critstack & bit) {
				active &= ~(stack1_54 << (3 * i)); // active must  critbf  
				if (minix[2] & mask)	AssignBf2(minix[2] & mask);
			}
		}
		active |= critbf;// be sure to keep it
	}

	inline void MinVmini(int imini, int vmini) {
		uint32_t shift = 3 * imini, bit9 = 1 << imini, cmin = 1,
			mask = 7 << shift, vminishift = vmini << shift;
		if (vmini == 8) {
			minix[0] |= bit9;// mini triplet
			critbf |= mask;
		}
		else {
			pairs27 |= vminishift;
			uint32_t cc = _popcnt32(vmini);
			minix[cc] |= bit9;
			if (cc > 1)critbf |= mask;
			else critbf |= (mask ^ vminishift);
			if (cc == 3) 	cmin = 2;
		}
		mincount += cmin;
		stackmin.bf.u16[imini % 3] += cmin;
	}
	inline int Stackx() {
		register uint32_t x = nb3;
		if (stackf.bf.u16[0] > x)return 1;
		if (stackf.bf.u16[1] > x)return 1;
		if (stackf.bf.u16[2] > x)return 1;
		if (stackf.bf.u16[0] == x)critstack |= 1;
		if (stackf.bf.u16[1] == x)critstack |= 2;
		if (stackf.bf.u16[2] == x)critstack |= 4;
		return 0;
	}
	inline void SetStackf(uint64_t bf) {
		stackf.bf.u64 = stackmin.bf.u64 + stackb12.bf.u64;
		nmiss = nb3 - mincount;
		if (!nmiss) critstack = 7;
	}
	inline int GetToAss() { return (nb3 - _popcnt32(assigned)); }
}scritb3;

int G17B::IsValid_myb12() {
	if (zh2b[1].IsValid(myb12)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 12) {
				cout << Char54out(ua54) << " seen add b12 size " << cc << endl;
				aigstop = 1; cout << " stop uab12<12" << endl; return 1;
			}
			if (cc < 22) {
				t54b12.AddA(ua54);
				t54b12.AddB(ua54);
				t54b12.AddC(ua54);
				t54b12.AddM(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
uint32_t G17B::IsValidB3(uint32_t bf) {
	if (zhou[0].CallCheckB3( bf)) {
		anduab3 = BIT_SET_27;
		stopexpandb3 = 0;
		//cout << "\t\t\t\t   add b3 nadd="<< zhgxn.nua << endl;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				cout << " bug no b3 IsValidB3 [7]" << p_cpt2g[7]
					<< "   [8]" << p_cpt2g[8] << endl;
				aigstop = 1;	return 1;
			}
			register uint32_t ua = w.bf.u32[2];
			t3_2[nt3_2++] = ua;
			anduab3 &= ua;
			register uint64_t U = w.bf.u64[0];
			uint64_t cc0 = _popcnt64(U);
			if (cc0 > 16)continue;
			U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
			if (cc <= 6 && cc>1) {
				if (cc > 3) 	myband3->Addguam(w, 1);				
				else {
					if (cc == 2) {
						int i81 = myband3->GetI81_2(w.bf.u32[2]);
						guah54.Add2(U, i81); guah54_2.Add2(U, i81);
						cb3.g2t.setBit(i81);
						p_cpt2g[20]++;
					}
					else {
						int i81 = myband3->GetI81_3(w.bf.u32[2]);
						guah54.Add3(U, i81); guah54_2.Add3(U, i81);
						cb3.g3t.setBit(i81);
						p_cpt2g[21]++;
					}
				}
			}
		}
		return 1;
	}
	else return 0;
}

void STD_B3::Go(CALLBAND3& cb3) {
	p_cpt2g[8]++;
	BF128 g2 = cb3.g2t & g.gsocket2,	g3 = cb3.g3t & g.gsocket3;
	register uint32_t  vmini, Mg2 = 0, Mg3 = 0;
	int x;
	while ((x = g2.getFirst128()) >= 0) {
		g2.clearBit(x);		Mg2 |= g.ua2bit27[x];	}
	while ((x = g3.getFirst128()) >= 0) {
		g3.clearBit(x);		Mg3 |= 1<<g.ua3_imini[x];	}
	scritb3.Init(cb3.ncl, cb3.cbs);
	for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
		if (!(vmini = Mg2 & 7))	if (Mg3 & 1)vmini = 8;
		if (vmini)scritb3.MinVmini(imini, vmini);
	}
	if ((int)scritb3.mincount >   scritb3.nb3) return;
	scritb3.SetStackf(cb3.bf12);// stackf; nmisss; critstack if nmiss=0
	if (scritb3.Stackx())return;// never > nb3 in stack
	scritb3.AssignCritical();// can assign 2 pairs in critical stacks
	if (!guam2done) BuildGuam2(g17b.myandall);

	// now get still valid uas band 3 not gua2s gua3s
	g17b.nt3 = 0;
	register uint64_t F = cb3.bf12;
	register uint32_t As = scritb3.assigned,Ac= scritb3.active;
	uint32_t t2[100], nt2 = 0, t3[100], nt3 = 0, t4[400], nt4 = 0, tm[400], ntm = 0;

	GUAM* mygu = (guam9done) ? tguam9 : tguam2;
	uint32_t ilim = (guam9done) ? ntguam9 : ntguam2;
	for (uint32_t i = 0; i < ilim; i++) {
		GUAM gm = mygu[i];
		if (gm.bf12 & F) continue;// assigned in bands 1+2
		register uint32_t U = gm.bf3;
		if (U & As) continue;// assigned in bands 3
		if (!(U &= Ac)) return; // dead branch
		int cc = _popcnt32(U);
		switch (cc) {
		case 1:g17b.t3[g17b.nt3++] = U; break;
		case 2:t2[nt2++] = U; break;
		case 3:t3[nt3++] = U; break;
		case 4:t4[nt4++] = U; break;
		default:tm[ntm++] = U; break;
		}		
	}
	for (uint32_t i = 0; i < nua; i++) {// now band 3 uas
		register uint32_t U = tua[i];
		if (U & As) continue;// assigned in bands 3
		if (!(U &= Ac)) return; // dead branch
		int cc = _popcnt32(U);
		switch (cc) {
		case 1:g17b.t3[g17b.nt3++] = U; break;
		case 2:t2[nt2++] = U; break;
		case 3:t3[nt3++] = U; break;
		case 4:t4[nt4++] = U; break;
		default:tm[ntm++] = U; break;
		}
	}
	// can now add all 2  if no redundancy
	for (uint32_t i = 0; i < nt2; i++) g17b.AddT3If(t2[i]);
	{
		register uint32_t V = scritb3.pairs27,i27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;
			g17b.AddT3If(g.pat2[i_27_to_81[i27]]);
		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				g17b.AddT3If (7 << (3 * i));
		}
	}
	// add now more clues 
	for (uint32_t i = 0; i < nt3; i++) g17b.AddT3If(t3[i]);
	for (uint32_t i = 0; i < nt4; i++) g17b.AddT3If(t4[i]);
	for (uint32_t i = 0; i < ntm; i++) g17b.AddT3If(tm[i]);
	memcpy(&genb12.grid0[54], band0, sizeof band0);// used in brute force
	g17b.myband3 = this;
	g17b.GoB3CleanOne();
}

void G17B::GoB3CleanOne() {// assign all singles in uas
	while (1) {
		int is1 = 0;
		for (int i = 0; i < (int)nt3; i++) {
			register uint32_t U = t3[i];
			if (_popcnt32(U) == 1)  is1 |= U;			
		}
		if (is1) {
			if (scritb3.AddAssign(is1)) return ;// conflict
			register uint32_t AC = scritb3.active,	F = scritb3.assigned;
			int n = nt3; nt3 = 0;
			for (int i = 0; i < n; i++) {
				register uint32_t U = t3[i];
				if (!(U & F)) {
					if (!(U &= AC)) return; // dead branch
					t3[nt3++] = U ;
				}
			}
		}
		else break;
	}
	memcpy(&grid0[54], myband3->band0, sizeof genb12.bands3[0].band0);
	if (_popcnt32(scritb3.assigned) == scritb3.nb3) {// finished 
		if ((nt3_2 = nt3)) return; // not valid if still uas
		GoB3End(0);
		return;
	}
	// sort t3 on size 

	switch (scritb3.nmiss) {
	case 0: GoB3Miss0(); return;
	case 1: GoB3Miss1(); return;
	default:GoB3MissMore(); return;
	}
}
void G17B::GoB3Miss0() {// no out field here assigned expected
	// store final  and check redundancy
	nt3_2 = 0;
	for (uint32_t i = 0; i < nt3; i++) {
		register uint32_t U = t3[i], nu = ~U;
		for (uint32_t j = 0; j < nt3_2; j++) {
			if (!(t3_2[j] & nu)) { U = 0; break; }
		}
		if (U)t3_2[nt3_2++] = U;
	}
	p_cpt2g[10]++;
	GoB3End(scritb3.GetToAss());
}
void G17B::GoB3Miss1() {// if outfield must have one out 
	int n2 = 0;
	{// store final  get outfield and check redundancy
		nt3_2 = 0;
		register uint32_t  critbf = scritb3.critbf,
			nmiss = scritb3.nmiss, n1 = 0;
		register uint32_t andout = scritb3.active & (~critbf);
		for (uint32_t i = 0; i < nt3; i++) {
			register uint32_t U = t3[i], nu = ~U;
			if (!(U & critbf)) { andout &= U; n2 = 1; continue; }
			for (uint32_t j = 0; j < nt3_2; j++) {
				if (!(t3_2[j] & nu)) { U = 0; break; }
			}
			if (U)t3_2[nt3_2++] = U;
		}
		if (n2) {
			if (!andout)return;// no possibility to add
			t3_2[nt3_2++] = andout; // dummy ua hitting all "out"
		}
	}
	if (!clean_valid_done) {// stop if not a valid band 1+2
		clean_valid_done = 1;
		if (IsValid_myb12()) { clean_valid_done = 2; return; }
	}
	p_cpt2g[11]++;
	if ((!n2) || (!scritb3.minix[2])) {// nothing to assign
		GoB3End(scritb3.GetToAss());
		return;
	}
	// now all critical can assign bf2
	//scritb3.critstack = 7;// not here in conflict with add one
	scritb3.AssignCritical();
	scritb3.active |= t3_2[nt3_2 - 1];// outfield lost to reenter 
	// reduce uas and look for fresh singles
	while (1) {
		int is1 = 0;
		register uint32_t AC = scritb3.active, F = scritb3.assigned;
		int n = nt3_2; nt3_2 = 0;
		for (int i = 0; i < n; i++) {
			register uint32_t U = t3_2[i];
			if (!(U & F)) {
				if (!(U &= AC)) return; // dead branch
				if (_popcnt32(U) == 1)  is1 |= U;
				t3_2[nt3_2++] = U;
			}
		}
		if (is1) {
			if (scritb3.AddAssign(is1)) return;// conflict
		}
		else break;
	}

	if (_popcnt32(scritb3.assigned) == scritb3.nb3) {
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (clean_valid_done == 2) return;
		if (IsValidB3(scritb3.assigned)) return;
		//_________________________   valid 18 
		Out17(scritb3.assigned);
		return;
	}
	GoB3End(scritb3.GetToAss());
}
void G17B::GoB3MissMore() {// no out field analysis here 
	// store final  and check redundancy
	p_cpt2g[12]++;
	nt3_2 = 0;
	for (uint32_t i = 0; i < nt3; i++) {
		register uint32_t U = t3[i], nu = ~U;
		for (uint32_t j = 0; j < nt3_2; j++) {
			if (!(t3_2[j] & nu)) { U = 0; break; }
		}
		if (U)t3_2[nt3_2++] = U;
	}
	GoB3End(scritb3.GetToAss());
}

void  G17B::GoB3End(int ntoass) {
	if (ntoass > 1) {GoB3Expand(ntoass); return;}
	if (!ntoass) {// can not have pending uas check final
		if (nt3_2) return;// should never be
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (clean_valid_done == 2) return;
		if (IsValidB3(scritb3.assigned)) return;
		Out17(scritb3.assigned);
		return;
	}
	if (ntoass != 1) {
		cout << "bug gob3end  not valid ntoass="<<ntoass 
			<< " [4] " << p_cpt2g[4] << " [7] " << p_cpt2g[7] << endl;
		aigstop = 1;		return;
	}
	// only one to ass must hit all uas 
	register int A = t3_2[0];
	for (int i = 1; i < (int)nt3_2; i++)
		if (!(A &= t3_2[i])) return;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValid_myb12()) { clean_valid_done = 2; return; }
	}
	register int cell;
	while (bitscanforward(cell, A)) {
		CRITB3 critb3 = scritb3;
		register uint32_t bit = 1 << cell;
		A ^= bit; //clear bit
		if (critb3.Addone(cell)) 	continue;
		if (IsValidB3(critb3.assigned)) {
			if (stopexpandb3) return;
			A &= anduab3;
		}
		else 	Out17(critb3.assigned);
	}


}

void G17B::GoB3Expand(int ntoass) {
	uint64_t limspot = ntoass - 1, limm1 = limspot - 1;
	struct SPB {
		CRITB3 critb3;
		uint32_t  possible_cells, indtw3;
	}spb[12];

	register SPB* s, * sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->critb3 = scritb3;
	s->indtw3 = 0;// initial first ua
	s->possible_cells = t3_2[0];
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)if (--s >= spb)goto next; else return;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->critb3.active ^= bit;
		sn = s + 1; *sn = *s;
		if (sn->critb3.Addone(cell)) 	goto next;		
	}
	if (ispot < limm1) {// here max 16 clues never valid b3
		register uint32_t F = sn->critb3.assigned;
		for (uint32_t i = s->indtw3 + 1; i < nt3_2; i++) {
			register uint32_t U = t3_2[i];
			if (!(U & F)) {
				U &= s->critb3.active;
				if (!U)goto next;//dead branch
				sn->possible_cells = U;
				sn->indtw3 = i;
				s++; // switch to next spot
				goto next;
			}
		}
		// no ua available must check( in 18 mode can not be valid)
		int rn = nt3_2; //first add if any
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			anduab3= t3_2[rn]& s->critb3.active;// first ua active cells
			//s->possible_cells &= anduab3;//  no this is a bug
			sn->possible_cells = anduab3;
			sn->indtw3 = rn;
			s++;	goto next;// switch to next spot
		}
		cout<< Char27out(sn->critb3.assigned) << "bug seen valid 16 or less ntoass= "
			<<ntoass<< endl;
		aigstop = 1;
		return;
	}


	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->critb3.active;
		register uint32_t F = sn->critb3.assigned;
		{// and still valid   uas
			for (uint32_t i = s->indtw3 + 1; i < nt3_2; i++) {
				register uint32_t U = t3_2[i];
				if (!(U & F)) { // not hit ua
					if (!(andw &= U))goto next;//dead branch
					aig = 0;
				}
			}
		}
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (aig) {// no ua could be a 17 valid in 18 mode 
			if (IsValidB3(F)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else { Out17(F);	goto next; }// this is a valid 17
		}
		register int cell;
		while (bitscanforward(cell, andw)) {
			CRITB3 critb3 = sn->critb3;
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else 	Out17(F | bit);			
		}
		goto next;
	}
	goto next;// never here
}

//=========brute force specific to this
int ZHOU::CallCheckB3( uint32_t bf, int nogo) {// 17 search mode
	zh_g.diag = nogo;
	memcpy(this, zhoustart, sizeof zhoustart);
	misc.SetAll_0();
	BF128 dca[9];
	memset(dca, 0, sizeof dca);
	int digitsbf = 0;
	{// cells in bands 1+2		
		int cell;
		register uint64_t U = g17b.myb12;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			int  digit = g17b.grid0[cell];
			int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
			digitsbf |= 1 << digit;
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = g17b.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zhgxn.nua = 0;
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	if(nogo)ImageCandidats();
	int ir = Full17Update();
	if (nogo)ImageCandidats();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);
	return zhgxn.nua;
}

int ZHOU::PartialInitSearch17(uint32_t* t, int n) {
	zh_g2.digitsbf = 0;
	memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = zh_g2.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		zh_g2.digitsbf |= 1 << digit;
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g2.Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | zh_g2.Digit_cell_Assigned[i];
	return 0;
}

int ZHOU::Apply17SingleOrEmptyCellsB12() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint64_t R1, R2, R3, R4;
	{
		register uint64_t* P = FD[0][0].bf.u64, M = *P;
		R1 = M;
		P += 4; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 4; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u64[0] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u64[0]; // these are new singles
	if (R1) {// singles to apply
		uint32_t xcell;
		while (bitscanforward64(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true singles
	if (!R2) {
		R3 &= ~R4;
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward64(zh_g2.xcell_to_guess, R2);
	return 0;
}
int ZHOU::Apply17SingleOrEmptyCellsB3() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint32_t R1, R2, R3, R4;
	{
		register uint32_t* P = &FD[0][0].bf.u32[2], M = *P;
		R1 = M;
		P += 8; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 8; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u32[2] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u32[2]; // these are new singles
	if (R1) {// singles to apply in band 3
		uint32_t xcell;
		while (bitscanforward(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = xcell + 54;
			xcell += 64;
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true pairs
	if (!R2) {
		R3 &= ~R4; // now true tripletss
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward(zh_g2.xcell_to_guess, R2);
	zh_g2.xcell_to_guess += 64;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		//if(diag &&cells_unsolved.bf.u32[2]==0)ImageCandidats();
		if (!Unsolved_Count()) return 2;
		if (cells_unsolved.bf.u32[2]) {// fill B3 first
			if (Apply17SingleOrEmptyCellsB3())	return 0; //  empty cell or conflict singles in cells
		}
		else {
			if ((!ISFALSEON))return 0;
			if (Apply17SingleOrEmptyCellsB12())	return 0;
		}
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Compute17Next(int index) {
	int ir = Full17Update();
	if (zh_g.diag) {
		cout << "index=" << index << endl;
		ImageCandidats();
	}
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128  wua;
			int* sol = g17b.grid0;
			wua.SetAll_0();
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {
				//if (zh_g.diag) {
					//cout << "zhgxn.nua=" << zhgxn.nua << endl;
					//wua.Print3(" ");
				//}
				int cc = _popcnt32(wua.bf.u32[2]);
				if ((!zhgxn.nua) || cc < 3)
					zhgxn.tua[zhgxn.nua++] = wua;
				if (cc < 3 || zhgxn.nua>5)	zh_g.go_back = 1;
			}
		}
		return;
	}
	Guess17(index);// continue the process
}
void ZHOU::Guess17(int index) {
	if (zh_g.go_back) return;
	int xcell = zh_g2.xcell_to_guess,
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU* mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1);
		if (zh_g.go_back) return;
	}
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (ISFALSEON && zhgxn.nua) continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU* mynext = (this + 1);
			*mynext = *this;
			if (cell >= 54)mynext->ISFALSEON++;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1);
			if (zh_g.go_back) return;
		}
	}
}

void G17B::Out17(uint32_t bfb3) {
	if (op.out_one) if (myband3->poutdone++) return;	
	char ws[128];
	strcpy(ws, empty_puzzle);
	{// cells in bands 1+2		
		int cell;
		register uint64_t U = myb12;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			ws[cell] = grid0[cell] + '1';
		}
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[i+54]= grid0[i+54] + '1';

	sprintf(&ws[81], ";%5d;%3d;%3d;%3d",(int)( genb12.nb12 / 64), genb12.i1t16, genb12.i2t16, t416_to_n6[myband3->i416]);     
	fout1 << ws << endl;
	a_17_found_here++;

}