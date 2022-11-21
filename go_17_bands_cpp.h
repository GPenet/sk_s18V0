

// standard first band (or unique band)

void STD_B416::Initstd() {
	strcpy(band, "12345678945");
	for (int i = 0; i < 11; i++) band0[i] = band[i] - '1';
	for (int i = 0; i < 27; i++)map[i] = i;// no mapping
	dband = 0;
}
void STD_B416::GetBandTable(int i416e) {
	i416 = i416e;
	strncpy(&band[11], t416[i416e], 16);
	band[27] = 0;
	for (int i = 11; i < 27; i++) band0[i] = band[i] - '1';
}
void STD_B416::SetGangster() {
	memset(gangster, 0, sizeof gangster);

	for (int ir = 0, i = 0; ir < 3; ir++)
		for (int ic = 0; ic < 9; ic++, i++)
			gangster[ic] |= 1 << band0[i];
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	memset(dpat, 0, sizeof dpat);
	memset(dpmi, 0, sizeof dpmi);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = band0[cell];
			dpmi[dig] |= Zhoucol << i;
			dpat[dig] |= 1 << cell;
			// compatibility to old code
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}
}
void STD_B416::InitC10(int i) { // band 1 mode enum
	GetBandTable(i); SetGangster();
	//zh1b_g.GetBand(band0, tua);// set zhone_i
}
void STD_B416::InitG12(int i) {// band 1 after initstd
	GetBandTable(i); SetGangster(); GetUAs();
	MorphUas(); // to add the count
}
void STD_B416::MorphUas() {
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i]&BIT_SET_27, ua = 0;
		register uint32_t cc;
		while (bitscanforward(cc, uao)) {
			uao ^= 1 << cc;
			ua |= 1 << map[cc];
		}
		tua[i] = ua | _popcnt32(ua) << 27;// store with count
	}

}
void STD_B416::InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
,int iband) {
	i416 = i16;
	dband = 27*iband;
	strncpy(band, ze, 27);
	band[27] = 0;
	for (int i = 0; i < 27; i++) band0[i] = band[i] - '1';
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
		map[vr0 + j] = vr + p.cols[j];
	}
	GetUAs();	MorphUas();	SetGangster();
}

void STD_B3::InitBand3(int i16, char* ze, BANDMINLEX::PERM& p) {
	InitBand2_3(i16, ze, p);
	//memset(&guas, 0, sizeof guas);
	memset(&g, 0, sizeof g);
	ntguam = poutdone=0;
	// setup minirows bit fields
	for (int i = 0; i < 9; i++) {
		minirows_bf[i] = 0;
		int* p = &band0[3 * i];
		for (int j = 0; j < 3; j++)
			minirows_bf[i] |= 1 << p[j];
	}
	//cout << band << " " << oct << minirows_bf[0]
	//	<< " " << minirows_bf[1] << dec << endl;

}

void STD_B416::PrintStatus() {
	cout << "band status i=" << i416 << "\tstart=" << dband << endl << "map ";
	for (int i = 0; i < 27; i++)cout << map[i] << " ";
	cout << endl;
	cout << band << endl << "gangster status" << endl;
	cout << "UAs table" << endl;
	for (uint32_t i = 0; i < nua; i++)
		cout << Char27out(tua[i]) << endl;
}


