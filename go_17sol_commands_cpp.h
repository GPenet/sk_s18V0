
//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };

const char * libs_c17_00_cpt2g[100] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 steps external loop ",//2
	"3 3 clues ",	
	"4 6 clues ",	
	"5 7p clues last ",
	"6 active 6 clues",
	"7 set b12 ",
	"8 go b3",
	"9 min too high",
	"10 miss0 ",
	"11 miss1 ",
	"12 missmore", 
	"13 missmore big","14 ","15","16 ","17 ","18 ","19 ",
	"20 addg2", 
	"21 addg3","22 ","23 ","24 ","25 ","26  ","27  ","28  ","29  ",	
	"30 nb3 min",
	"31 nb3 max",
	"32 nb3 tot  ","33 ","34 ",
	"35 uas min",
	"36 uas max",
	"37 uas tot","38 ", "39 ",
	"40 ng2 min",
	"41 n2 max",
	"42 ng2 tot","43 ","44 ",
	"45 ng3min",	
	"46 ng3 max",
	"47 ng3 tot","48 ","49 ",
	"50 ",
	"51 full","52 below",
	"53 nfull<20",
	"54 size 7",
	"55 size 8",	
	"56 size 9",
	"57 size 10",
	"58 size 11",
	"59 ",
	"60 max 12 clues",
	"61 max 11 clues",
	"62 max 10 clues","63 ","64 ","65 ",	"66 ","67 ","68 ","69 ",
	"70 count entry go 7_12",
	"71 count 70 not empty",
	"72 max 7_9 ",
	"73 count 12 direct",
	"74 max 12 direct",
	"75 ",	"76 ","77 ","78 ","79 ",
	"80 ","81","82 ","83 ","84 ","85 ",	"86 ","87 ","88 ","89 ",
	"90 ","91","92 ","93 ","94 ","95 ",	"96 ","97 ","98 ","99 ",

};
void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17 clues  " << endl;
	op.SetUp(0);
	int it16_start = sgo.vx[0];
	g17b.aigstop=0;
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt1g, 0, sizeof p_cpt1g);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(op.b1);
	cout << "print final stats" << endl;
	for (int i = 0; i < 100; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "exit final stats" << endl;
}
//========================= known s18 
void Go_c17_10( ) {
	//op.SetUp(10,0);
	cout << "Go_10() search 17/18 using a file having known 17 656 " << endl;
	op.SetUp(10,1);// setup with known
	cout << "back setup " << endl;

	zh_g.modevalid = 1;
	zh_g2.grid0 = genb12.grid0;
	zh_g2.zsol = zh_g2.stdfirstsol;
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	char* ze2 = &ze[82];
	//if (op.t18) return;
	//if (op.p1) return;
	//cout << "Go_10() search 17/18 using a file having known 17 656 " << endl;

	//op.ton = 1;

	while (finput.GetLigne()) {
		if(strlen(ze)<160) continue;// skip blank lines
		npuz++;
		g17b.npuz = npuz;
		g17b.a_17_found_here = 0;
		g17b.aigstop= 0;
		if (npuz <= op.skip) continue;
		if (npuz > op.last) break;
		CBS cbs;
		g17b.p17diag.SetAll_0();
		int  nclues = 0;
		memset(&cbs, 0, sizeof cbs);
		for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
			g17b.p17diag.Set_c(i);
			cbs.Add(i);
			nclues++;
		}

		if (op.ton)cout << ze << " to process  n=" << npuz
			<< " bands " << cbs.b[0] << cbs.b[1] << cbs.b[2]
			<< " stacks " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
		if (op.t18 && nclues == 17) {
			if (op.ton)cout << " see later 17 clues for a 18" << endl;
			continue;
		}
		if (op.p2 && (cbs.s[0] > 6 || cbs.s[1] > 6 || cbs.s[2] > 6)) {
			if (op.ton)cout << " pass2 stack > 6 clues" << endl;
				continue;
		}
		if (op.t18) {
			if (op.p2 && (cbs.b[0] > 6 || cbs.b[1] > 6 || cbs.b[2] > 6)) {
				if (op.ton)cout << "t18  pass2 not 666 666" << endl;
					continue;
			}
			if (!op.p2) {
				if (cbs.b[2] < 7) {
					if (op.ton)cout << "t18  pass1 not >=7 in b3" << endl;
						continue;
				}
				if (cbs.b[2] < cbs.s[0] || cbs.b[2] < cbs.s[1] || cbs.b[2] < cbs.s[2] ) {
					if (op.ton)cout << "t18  stack too high" << endl;
						continue;
				}
			}
		}
		if ((!op.t18) && op.p2) {
			op.p2b = 0;
			//================================ to avoid the 665 case
			if (cbs.b[2] == 5) {// change band3 <-> band2
				op.p2b = 1;
				for (int i = 0; i < 27; i++) {
					char temp = ze[i + 27];	ze[i + 27] = ze[i + 54];	ze[i + 54] = temp;
					temp = ze[i + 109];	ze[i + 109] = ze[i + 136];	ze[i + 136] = temp;
				}
				uint32_t w = g17b.p17diag.bf.u32[1];
				g17b.p17diag.bf.u32[1]= g17b.p17diag.bf.u32[2];
				g17b.p17diag.bf.u32[2]=w;
				
			}


		}

		cout << "\n\nto process  n="<<dec << npuz <<" debug="<< op.ton << endl;
		if (op.ton)		cout << ze << " to process  n="  << npuz 
			<< " bands "<< cbs.b[0] << cbs.b[1] << cbs.b[2] 
			<< " stacks " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;


		// =======================morph entry 
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		myband1.InitBand2_3(ib1, ze, perm_ret, 0);
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		myband2.InitBand2_3(ib2, &ze[27], perm_ret, 1);
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		genb12.bands3[0].InitBand3(ib3, &ze[54], perm_ret);
		genb12.nband3 = 1;

		ze[81] = 0;
		if (op.ton)
			cout << Char2Xout(g17b.p17diag.bf.u64[0]) << " b12 pattern for the 17" << endl;
		register uint64_t U = g17b.p17diag.bf.u64[0];
		g17b.pk54= (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);
		genb12.ValidInitGang();
		g17b.npuz = npuz;
		if(sgo.bfx[2] & 1)g17b.StartKnown();
		else	g17b.Start();
		//g17b.a_17_found_here = 1;
		if (!g17b.a_17_found_here) {
			cout << "puz="<<npuz << " failed to find the searched 17" << endl;
			cerr << "puz=" << npuz << " failed to find the searched 17" << endl;
			break;
		}
	}
	cout << "print final stats" << endl;
	for (int i = 0; i < 100; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}

void Go_c17_12() {// check diagonal status in a 665
	// search 17 using a file having known  as entry and one 17 given 6 6 5
	char * ze = finput.ze;
	int  zs0 [81],zs0d[81],
		npuz = 0;
	cout << "Go_c17_12() band analysis option=" <<sgo.vx[2] << endl;
	while (finput.GetLigne()) {
		if (strlen(ze) < 81)continue;
		// =======================morph entry to have min n6 count in first
		for (int i = 0; i < 81; i++) {
			zs0[i] = ze[i] - '1';
			zs0d[i] = ze[C_transpose_d[i]] - '1';
		}
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		int ib1 = perm_ret.i416, ib1a = t416_to_n6[ib1];
		bandminlex.Getmin(&zs0[27], &perm_ret);
		int ib2 = perm_ret.i416, ib2a = t416_to_n6[ib2];
		bandminlex.Getmin(&zs0[54], &perm_ret);
		int ib3 = perm_ret.i416, ib3a = t416_to_n6[ib3];
		bandminlex.Getmin(zs0d, &perm_ret);
		int ib1d = perm_ret.i416, ib1ad = t416_to_n6[ib1d];
		bandminlex.Getmin(&zs0d[27], &perm_ret);
		int ib2d = perm_ret.i416, ib2ad = t416_to_n6[ib2d];
		bandminlex.Getmin(&zs0d[54], &perm_ret);
		int ib3d = perm_ret.i416, ib3ad = t416_to_n6[ib3d];

		if (0) {
			cout << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
				<< ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< endl;
		}
		switch (sgo.vx[2]) {
		case 1:// extract 18 (8 in [87]
			if (strlen(ze) < 88)break;
			if (ze[87] == '8')
				fout1 << ze << endl;
			break;
		
		case 2:// extract 18 (8 in [87]) plus info
			if (strlen(ze) < 88)break;
			if (ze[87] == '8')
				fout1 << ze << ";\t" << ib1 << ";" << ib2 << ";" << ib3
				<< ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1d << ";" << ib2d << ";" << ib3d
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< endl;
			break;
		case 3: {// sol+18 to analyze get pass2
			if (strlen(ze) < 162)break;
			if (!sgo.vx[3])sgo.vx[3] = 82;
			char* ze2 = &ze[sgo.vx[3]];
			CBS cbs;
			int  nclues = 0;
			memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
				cbs.Add(i);
			}
			if(cbs.b[0] ==6 && cbs.b[1] ==6 &&  cbs.s[0]==6 && cbs.s[1]==6)
			fout1 << ze2 << ";\t"  << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< "; bands; " << cbs.b[0] << cbs.b[1] << cbs.b[2]
				<< "; stacks; " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;

		}
		case 4: {// sol+18 to analyze get pass1
			if (strlen(ze) < 162)break;
			if (!sgo.vx[3])sgo.vx[3] = 82;
			char* ze2 = &ze[sgo.vx[3]];
			CBS cbs;
			int  nclues = 0;
			memset(&cbs, 0, sizeof cbs);
			for (int i = 0; i < 81; i++) if (ze2[i] != '.') {
				cbs.Add(i);
			}

			if (cbs.b[0] == 6 && cbs.b[1] == 6 && cbs.s[0] == 6 && cbs.s[1] == 6)
				break;// pass2
			uint16_t x=cbs.s[0];
			for (int i = 1; i < 3; i++) if (x < cbs.s[i])x = cbs.s[i];
			if (cbs.b[2] < cbs.b[0]) break;
			if (cbs.b[2] < cbs.b[1]) break;
			if (cbs.b[2] < x) break;
			fout1 << ze2 << ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< "; bands; " << cbs.b[0] << cbs.b[1] << cbs.b[2]
				<< "; stacks; " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;
			if(cbs.b[0] + cbs.b[1] ==7)
				cout << ze2 << ";\t" << ib1a << ";" << ib2a << ";" << ib3a
				<< ";\t" << ib1ad << ";" << ib2ad << ";" << ib3ad
				<< "; bands; " << cbs.b[0] << cbs.b[1] << cbs.b[2]
				<< "; stacks; " << cbs.s[0] << cbs.s[1] << cbs.s[2] << endl;

		}
		}
	}

}
