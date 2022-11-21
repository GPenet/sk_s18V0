/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
//====================================
#include "sk_t.h"  // tab0  
#include "Zhn.h"
#include "Zhtables_cpp.h" // also describes the pattern
//#define DIAG
ZH_GLOBAL2 zh_g2;
 /*ZH_1D class to solve one digit 3 bands
all valid solutions  are stored in table
the table is supplied by the caller
this is somehow a private class for the brute force
and this is part of the critical code in the brute force 
except for easiest puzzles immediatly solved
struct ZH_1D_GLOBAL {
	BF128 *tsolw,t3; // belong to the caller
	int nsolw;
	ZH_1D_GLOBAL() { t3.SetAll_1(); t3.bf.u32[3] = 0; }
	inline void Add(BF128 & s) {
		*tsolw++ = s & t3; // clear last 32 bits
		nsolw++;
	}
	int Go(BF128 & fde, BF128 *tsol);
};*/

struct ZH_1D {
	BF128 FD;// last 32 bits contains unsolved rows
	//int GetSols(int ru);
	//int GetAllSols(BF128 & fde, int ru);
	inline void Assign(int cell) {
		FD &= AssignMask_Digit[cell];
	}
	void Guess();
	void GuessGo(int cell);
	int Update();
};

ZH_1D_GLOBAL zh1d_g;
ZH_1D zh1d[10];
// _______here to try to optimize the cache the ZH_1D code
int ZH_1D_GLOBAL::Go(BF128 & fde, BF128 *tsol, int limit) {
	if (0)cout << "go 1d limit=" << limit << endl;
	lim = limit;
	tsolw = tsol;
	nsolw = 0;
	zh1d[0].FD = fde;
	zh1d[0].Guess();
	return nsolw;
}
void ZH_1D::Guess() {// not yet solved
	if (zh1d_g.nsolw > zh1d_g.lim) return;
	uint32_t rmask[3] = { 0777,0777000,0777000000 };
	uint32_t ru = FD.bf.u32[3], bu = TblShrinkMask[ru],
		minband = 30, minbandindex;
	if (0) {
		cout << "index=" << this - zh1d<< "nperms="<< zh1d_g.nsolw << endl;
		cout << Char9out(ru) << " ru" << endl;
		cout << Char9out(bu) << " bu" << endl;
	}
	// use the band with the lowed count
	for (int ib = 0, bit = 1; ib < 3; ib++, bit <<= 1) {
		if (bu&bit) {// not solved band
			uint32_t cc = _popcnt32(FD.bf.u32[ib]);
			if (cc <= 5) { minbandindex = ib; break; }
			if (cc < minband) {
				minband = cc;
				minbandindex = ib;
			}
		}
	}

	if (0) 		cout << "guess minbandindex=" << minbandindex << endl;
	// now select the unsolved  row with the lowest count
	uint32_t band = FD.bf.u32[minbandindex],
		bru = (ru >> (3 * minbandindex)) & 7,
		dcell = 27 * minbandindex, count = 10, myrow, myirow;
	for (int irow = 0, bit = 1; irow < 3; irow++, bit <<= 1) {
		if (bru&bit) {// not solved row in band
			uint32_t row = band & rmask[irow];
			uint32_t cc = _popcnt32(row);
			if (cc < count) {
				count = cc;
				myirow = irow;
				myrow = row;
				if (cc == 2) break;
			}
		}
	}
	if (0) 		cout << "guess rowindex=" << myirow << endl;
	// try each candidate in the selected row
	FD.bf.u32[3] ^= 1 << (3 * minbandindex + myirow);// set row solved
	uint32_t cell;
	while (bitscanforward(cell, myrow)) {
		myrow ^= 1 << cell;
		cell += dcell;
		(this + 1)->GuessGo(cell);
		if (zh1d_g.nsolw > zh1d_g.lim) return;
	}
}
void ZH_1D::GuessGo(int cell) {
	*this = *(this - 1);
	FD &= AssignMask_Digit[cell];// unsolved row already updated
	//int ru = FD.bf.u32[3];
	if (Update()) return; // locked
	//if (!(FD.bf.u32[3] = ru)) zh1d_g.Add(FD);
	if (!FD.bf.u32[3] ) zh1d_g.Add(FD);
	else Guess();
}
int ZH_1D::Update() {
	//char ws[82];
	//cout << FD.String3X(ws) << " debut update" << endl;
	register int Shrink, S, A, B, C, D, ru = FD.bf.u32[3],
		rub = TblShrinkMask[ru];// rub 3 bits 0 to 111
	if (0)cout << "update1d ru=0" << oct << ru
		<< " rub" << rub << dec << endl;

	switch (rub) {
	case 1: {// band 1 alone not solved
		A = FD.bf.u32[0];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[0] = A;
		ru = (ru & 0770) | S;
		break;// no reason to loop here
	}
	case 2: {// band 2 alone not solved
		A = FD.bf.u32[1];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[1] = A;
		ru = (ru & 0707) | (S << 3);
		break;// no reason to loop here
	}
	case 4: {// band 3 alone not solved
		A = FD.bf.u32[2];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[2] = A;
		ru = (ru & 077) | (S << 6);
		break;// no reason to loop here
	}
	case 3:while (1) {// band 1;2 not solved
		A = FD.bf.u32[0];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[1] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[0] = B = A;
		ru = (ru & 0770) | S;
		A = FD.bf.u32[1];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[0] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[1] = A;
		ru = (ru & 0707) | (S << 3);
		if (B == FD.bf.u32[0])return 0;;// no reason to loop  
	}
	case 5: while (1) {// band 1;3 not solved
		A = FD.bf.u32[0];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[2] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[0] = B = A;
		ru = (ru & 0770) | S;
		A = FD.bf.u32[2];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[0] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[2] = A;
		ru = (ru & 077) | (S << 6);
		if (B == FD.bf.u32[0])	return 0;// no reason to loop  
	}

	case 6: while (1) {// band 2;3 not solved
		A = FD.bf.u32[1];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[2] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[1] = B = A;
		ru = (ru & 0707) | (S << 3);
		A = FD.bf.u32[2];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[1] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[2] = A;
		ru = (ru & 077) | (S << 6);
		if (B == FD.bf.u32[1])	return 0;// no reason to loop  
	}
	case 7: {// no band  solved
		//cout << "loop1" << endl;
	loop1:
		A = FD.bf.u32[0];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		B = FD.bf.u32[2];
		FD.bf.u32[1] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		B = FD.bf.u32[1];
		FD.bf.u32[2] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[0] = C = A;
		ru = (ru & 0770) | S;
	loop2:
		//cout << "loop2" << endl;
		//cout << FD.String3X(ws) << " après loop1" << endl;
		A = FD.bf.u32[1];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		B = FD.bf.u32[2];
		FD.bf.u32[0] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		B = FD.bf.u32[0];
		FD.bf.u32[2] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[1] = D = A;
		ru = (ru & 0707) | (S << 3);
		//cout << FD.String3X(ws) << " après loop2" << endl;
		A =  FD.bf.u32[2];
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 1;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		B = FD.bf.u32[1];
		FD.bf.u32[0] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		B = FD.bf.u32[0];
		FD.bf.u32[1] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		FD.bf.u32[2] = A;
		ru = (ru & 077) | (S << 6);
		if (C != FD.bf.u32[0])goto loop1;
		if (D != FD.bf.u32[1])goto loop2;
	}
	}// end switch
	FD.bf.u32[3] = ru;
	//cout << FD.String3X(ws) << " sortie update" << endl;
	//if (1)cout << "exit update1d ru=0" << oct << ru	 << dec << endl;
	return 0;
}
//============================= ZH_GLOBAL code and workinfg areas

ZH_GLOBAL zh_g;
ZHOU  zhou_ip,//
zhou_solve;// basis to solve a puzzle using elimination logic
ZHOU zhou[50]; // must host main brute force plus minimality analysis and recursive generation
 ZH_GLOBAL::ZH_GLOBAL(){
	diag = 0;
	modevalid = 0;
	modeguess = 1;
	zh_g2.zsol =  0; // no solution unless required buy the user
	// init Tblgame
	for (int i = 0; i < 3; i++)init_3x.bf.u32[i] = BIT_SET_27;
	init_3x.bf.u32[3] = 0;
	init_digit = init_3x;
	init_digit.bf.u32[3] = 0777;//all rows unsolved
}

int ZH_GLOBAL::Go_InitSudoku(char * ze){
	ze[81] = 0;
	strcpy(zh_g2.puz, ze);
	zh_g2.ngiven = 0;
	int digs = 0;
    for (int i = 0; i < 81; i++){
		register int c = ze[i];
		if (c<'1' || c>'9') continue;
		c -= '1';
		digs |= 1 << c;
		zh_g2.tgiven[zh_g2.ngiven++].u16 =(uint16_t)( i | (c << 8));
	}
	if (_popcnt32(digs) < 8) return 1; // don't accept less than 8 digits given
	return zhou[0].InitSudoku(zh_g2.tgiven, zh_g2.ngiven);
}
void ZH_GLOBAL::ValidPuzzle(ZHOU * z){
	nsol++;
	if (!modevalid) {// standard process
		if (zh_g2.zsol && (nsol==1)) {// store the first solution
			z->SetKnown(zh_g2.zsol);
			for (int i = 0; i < 9; i++)zh_g2.digit_sol[i] = z->FD[i][0];
		}
	}
	else if (modevalid==1) {// 17 search stop at first solution grid
		z->SetKnown(zh_g2.zsol);
		nsol++;
		go_back = 1;
		return;
	}
	if (nsol > lim)go_back = 1;
}

//============================= zhou code

int ZHOU::Isvalid() { // usually after init 2 steps
	zh_g.Init(1);// maxsols=1
	ComputeNext();  return zh_g.nsol;
}
int ZHOU::CheckValidityQuick(char *puzzle) {
	zh_g.Init(1);// maxsols=1
	if (zh_g.Go_InitSudoku(puzzle)) return 0;
	int ir = FullUpdate();
	if (!ir) return 0;// not valid
	if (ir == 2) {
		if(zh_g2.zsol)zh_g.ValidPuzzle(this);
		return 1;// solved can not be multiple
	}
	Guess();
	return zh_g.nsol;
}
int ZHOU::FullUpdate() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
int ZHOU::ApplySingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here only singles and empty cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0], R3; 	R1 |= FD[1][0];
	R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs
		zh_g.pairs = R2 - R3;
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
void ZHOU::GuessGo(int dig, BF128& s) {
	*this = *(this - 1);
	BF128 assigned_cells = s & cells_unsolved;
	FD[0][0] -= assigned_cells;	FD[1][0] -= assigned_cells;
	FD[2][0] -= assigned_cells;	FD[3][0] -= assigned_cells;
	FD[4][0] -= assigned_cells;	FD[5][0] -= assigned_cells;
	FD[6][0] -= assigned_cells;	FD[7][0] -= assigned_cells;
	FD[8][0] -= assigned_cells;	cells_unsolved -= assigned_cells;
	FD[dig][0] = s;// restore digit pm
	if (0) {
		cout << "GuessV3Go digit " << dig + 1 << endl;
		//ImageCandidats();
	}
	ComputeNext();
}
void ZHOU::Guess() { 
	if (zh_g.go_back) return;
	if (cells_unsolved.isEmpty()) {
		//if (zh_g.diag)cout << endl << "valid guess v3\n" << endl;
		zh_g.ValidPuzzle(this);
		return;
	}
	if (zh_g.pairs.isNotEmpty()) {	GuessInCell();	return;	}
	if (GuessHiddenBivalue()) return;
	// no pair, no bi valuesolve a full digit 
	//cout << "guess full digit" << endl;
	GuessFullDigit();
}
void ZHOU::GuessInCell() {
	{ // select band with more unsolved cells 
		uint32_t nfreecells = 0, nw;
		if (zh_g.pairs.bf.u32[0]) {
			nfreecells = _popcnt32(cells_unsolved.bf.u32[0]);
		}
		if (zh_g.pairs.bf.u32[1]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[1]);
				if (nw > nfreecells) {
					nfreecells = nw;
					zh_g.pairs.bf.u32[0] = 0;
				}
			}
			else	nfreecells = _popcnt32(cells_unsolved.bf.u32[1]);
		}
		if (zh_g.pairs.bf.u32[2]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[2]);
				if (nw > nfreecells) {
					nfreecells = nw;
					zh_g.pairs.bf.u32[0] = 0;
					zh_g.pairs.bf.u32[1] = 0;
				}
			}
		}
	}
	int xcell = zh_g.pairs.getFirst128(),
		cell = From_128_To_81[xcell], tdig[2], ndig = 0;
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig][0].On(xcell))tdig[ndig++] = idig;
	ZHOU * mynext = this + 1; // start next guess
	*mynext=*this;
	mynext->SetaCom(tdig[0], cell, xcell);
	mynext->ComputeNext();
	if (zh_g.go_back) return;
	SetaCom(tdig[1], cell, xcell);
	ComputeNext();
}
int ZHOU::GuessHiddenBivalue() {// look a hidden pair in row or box
	uint32_t hidden;
	int idig;
	int dcell, dxcell;
	for (idig = 0; idig < 9; idig++) {
		BF128 & fd = FD[idig][0];
		register int Rows;
		if (!(Rows = fd.bf.u32[3])) continue;
		if (Rows & 7) {//try row band1
			dcell = dxcell = 0;
			register int  band = fd.bf.u32[0];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 070) {// try row band2
			dcell = 27; dxcell = 32;
			register int  band = fd.bf.u32[1];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 0700) {//try row band3
			dcell = 54; dxcell = 64;
			register int  band = fd.bf.u32[2];
			hidden = band & 0777;			if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0777000000;		if (_popcnt32(hidden) == 2)goto exitok;
		}
	}
	// no bi value row, try bi value box
	for (idig = 0; idig < 9; idig++) {// priority to high digits last done
		BF128 & fd = FD[idig][0];
		register int Rows;
		if (!(Rows = fd.bf.u32[3])) continue;
		if (Rows & 7) {//try bow band1
			dcell = dxcell = 0;
			register int  band = fd.bf.u32[0];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 070) {// try box band2
			dcell = 27; dxcell = 32;
			register int  band = fd.bf.u32[1];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
		if (Rows & 0700) {//try bow band3
			dcell = 54; dxcell = 64;
			register int  band = fd.bf.u32[2];
			hidden = band & 07007007;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 070070070;		if (_popcnt32(hidden) == 2)goto exitok;
			hidden = band & 0700700700;		if (_popcnt32(hidden) == 2)goto exitok;
		}
	}
	return 0;
exitok:
	uint32_t res;
	bitscanforward(res, hidden);
	ZHOU * mynext = this + 1; // start next guess
	*mynext=*this;
	mynext->SetaCom(idig, res + dcell, res + dxcell);
	mynext->ComputeNext();
	if (zh_g.go_back) return 1;
	bitscanreverse(res, hidden);
	SetaCom(idig, res + dcell, res + dxcell);
	ComputeNext();
	return 1;
}
void  ZHOU::GuessFullDigit() {
	// select the digit not fully solved lowest count of candidates
	int mincount = 90, digmin;
	for (int idig = 0; idig < 9; idig++) {
		int cc = (FD[idig][0] & zh1d_g.t3).Count();
		if (cc < 10)continue;
		if (cc < mincount) {
			mincount = cc; digmin = idig;
		}
	}
	if (zh_g.diag) {
		cout << " solve full digit=" << digmin + 1 << endl;
		ImageCandidats();
		Debug(1);
	}
	// collect perms valid for this digit
	BF128  tres[300];// to store possible perms 
	int nperms = zh1d_g.Go(FD[digmin][0], tres);
	if (zh_g.diag)cout << " nperms=" << nperms << endl;
	if (zh_g.diag &&nperms > 100) {
		cout <<"more than 100 perms " << nperms << "  digit=" << digmin + 1 << endl;
		ImageCandidats();
		Debug(1);
	}
	for (int i = 0; i < nperms; i++) {
		ZHOU * mynext = this + 1; // start next guess
		(this + 1)->GuessGo(digmin, tres[i]);
		if (zh_g.go_back) return;
	}
}
void ZHOU::Assign(int digit, int cell, int xcell) {
	FD[digit][0] &= AssignMask_Digit[cell];
	cells_unsolved.Clear(xcell);
}

void ZHOU::Setcell(int cell){ 
	SetaCom(zh_g2.zerobased_sol[cell], cell, C_To128[cell]);
}
void ZHOU::SetaCom(int digit, int cell, int xcell) { // single in cell
	BF128 *Fd = FD[digit];
	*Fd &= AssignMask_Digit[cell];
	cells_unsolved.Clear(xcell);
	BF128 * RF = FD[8];
	for (; RF >= FD[0]; RF -= 2)RF->Clear(xcell);
	Fd->setBit(xcell); // restore bit for digit assigned
}

#define UPD_AND_CL(W,X)	last_assigned = W+1;cur_assigned = 0;\
B = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[X] &= B;\
FD[0][0].bf.u32[X] &= B; FD[1][0].bf.u32[X] &= B; FD[2][0].bf.u32[X] &= B; \
FD[3][0].bf.u32[X] &= B; FD[4][0].bf.u32[X] &= B; FD[5][0].bf.u32[X] &= B; \
FD[6][0].bf.u32[X] &= B; FD[7][0].bf.u32[X] &= B; FD[8][0].bf.u32[X] &= B; }\
FD[W][0].bf.u32[X] = FD[W][1].bf.u32[X] = A;

#define UPD_012(W,X,Y,Z)A = FD[W][0].bf.u32[X];\
if (A != FD[W][1].bf.u32[X]){\
Shrink = (TblShrinkMask[A & 0x1FF]\
| TblShrinkMask[(A >> 9) & 0x1FF] << 3\
| TblShrinkMask[(A >> 18)] << 6);\
if ((A &= TblComplexMask[Shrink]) == 0)  return 0;\
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);\
B = FD[W][0].bf.u32[Y];\
FD[W][0].bf.u32[Z] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];\
B = FD[W][0].bf.u32[Z];\
FD[W][0].bf.u32[Y] &= TblMaskSingle[S] & TblMaskDouble[S | ((B | (B >> 9) | (B >> 18)) & 0x1FF)];\
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];

#define UPD_ONE_DIGIT(W) if (cur_assigned > W+1)goto exit_digits;\
	if (FD[W][0].bf.u64[0] != FD[W][1].bf.u64[0]\
|| FD[W][0].bf.u32[2] != FD[W][1].bf.u32[2]){\
r_free = FD[W][0].bf.u32[3];\
UPD_012(W, 0, 1, 2)	if ((r_free & 7) != S){\
r_free &= 0770 | S;	UPD_AND_CL(W, 0)}\
UPD_012(W, 1, 0, 2)	if (((r_free >> 3) & 7) != S){\
r_free &= 0707 | (S << 3);	UPD_AND_CL(W, 1)}\
UPD_012(W, 2, 0, 1)	if (((r_free >> 6) & 7) != S){\
r_free &= 077 | (S << 6);	UPD_AND_CL(W, 2)}\
FD[W][0].bf.u32[3] = r_free;}

int ZHOU::Update(){
	//if (zh_g.diag)cout << "Update index=" << index << endl;
	register uint32_t Shrink = 1, r_free, B, A, S, last_assigned = 0,cur_assigned;
loop_upd:
	//zh_g.cpt[3]++;
	cur_assigned = last_assigned; last_assigned = 0;
	if (FD[8][0].bf.u32[3]) { UPD_ONE_DIGIT(8) }
	if (FD[7][0].bf.u32[3]) { UPD_ONE_DIGIT(7) }
	if (FD[6][0].bf.u32[3]) { UPD_ONE_DIGIT(6) }
	if (FD[5][0].bf.u32[3]) { UPD_ONE_DIGIT(5) }
	if (FD[4][0].bf.u32[3]) { UPD_ONE_DIGIT(4) }
	if (FD[3][0].bf.u32[3]) { UPD_ONE_DIGIT(3) }
	if (FD[2][0].bf.u32[3]) { UPD_ONE_DIGIT(2) }
	if (FD[1][0].bf.u32[3]) { UPD_ONE_DIGIT(1) }
	if (FD[0][0].bf.u32[3]) { UPD_ONE_DIGIT(0) }
exit_digits:
	//if (zh_g.diag>1){ cout << "avant test loop last_assigned=" << last_assigned << endl;
	//Debug(1); ImageCandidats();
	//}
	if (last_assigned) goto loop_upd;// nothing to do in the next cycle
	return 1;
}

int ZHOU::InitSudoku(GINT16 * t, int n){ 
	BF128 Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | Digit_cell_Assigned[i];
	return 0;
}
int ZHOU::InitSudoku(char * zpuz){
	memmove(zh_g2.puz, zpuz, 82);
	return InitSudoku(zh_g2.tgiven, zh_g2.ngiven);
}
char * ZHOU::SetKnown(char * zs){
	strcpy(zs, empty_puzzle);
	BF128  fd;
	for (int digit = 0; digit < 9; digit++){
		fd = FD[digit][0]; // unsolved or partially solved, still in zhou
		if (fd.bf.u32[3] == 0777) continue;
		for (int ib = 0; ib < 3; ib++){// 3 blocs per digit
			int arows = (fd.bf.u32[3] >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = fd.bf.u32[ib];
			for (int j = 0; j<3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				uint32_t  irow;
				bitscanforward(irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = digit + '1';
			}
		}
	}
	return zs;
}

int ZHOU::IsMinimale(GINT16 * to, int no){// assumed checked valid before
	GINT16 td[81], cx;
	int nd;
	for (int i = 0; i < no; i++){
		nd = 0;
		for (int i2 = 0; i2 < no; i2++){
			if (i2 - i)td[nd++] = to[i2];
			else cx = to[i2];
		}
		zh_g.Init(0);
		InitSudoku(td, nd); // always valid
		ClearCandidate_c(cx.u8[1], cx.u8[0]);// clean the valid known
		//cout << "check minimale step" << endl;
		//ImageCandidats();
		ComputeNext();
		if (!zh_g.nsol) return 0; //could check also not valid
	}
	return 1;
}

#include "Zhn_doc_debug_cpp.h"
/*Zhn Control of flow
-b9- use of hidden pairs triplets zhg_modeguess
 abc  a use pairs  b use triplets c use pairs not hidden pair
 -v7- max index for hidden pairs ...
*/


int ZHOU::PartialInitSudoku(GINT16 * t, int n){// if morph, done before
	memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g2.Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | zh_g2.Digit_cell_Assigned[i];
	return 0;
}
int ZHOU::EndInitSudoku(GINT16 * t, int n){// if morph, done before
	BF128 Digit_cell_Assigned[9];
	memmove(Digit_cell_Assigned, zh_g2.Digit_cell_Assigned, 9 * 16);
	*this = zhou_ip;
	for (int ic = 0; ic < n; ic++)   {
		int digit = t[ic].u8[1], cell = t[ic].u8[0], xcell = C_To128[cell];
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved; w.bf.u32[3] = ~0;
	for (int i = 0; i<9; i++)  FD[i][0] &= w | Digit_cell_Assigned[i];
	return 0;
}


//====================== code for UA gen2 digits

#define UPD_AND_CL2(W,X)	last_assigned = W+1;cur_assigned = 0;\
B = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[X] &= B;\
FD[0][0].bf.u32[X] &= B; FD[1][0].bf.u32[X] &= B; } \
FD[W][0].bf.u32[X] = FD[W][1].bf.u32[X] = A;

#define UPD_ONE_DIGIT2(W) if (cur_assigned > W+1)goto exit_digits;\
	if (FD[W][0].bf.u64[0] != FD[W][1].bf.u64[0]\
|| FD[W][0].bf.u32[2] != FD[W][1].bf.u32[2]){\
r_free = FD[W][0].bf.u32[3];\
UPD_012(W, 0, 1, 2)	if ((r_free & 7) != S){\
r_free &= 0770 | S;	UPD_AND_CL2(W, 0)}\
UPD_012(W, 1, 0, 2)	if (((r_free >> 3) & 7) != S){\
r_free &= 0707 | (S << 3);	UPD_AND_CL2(W, 1)}\
UPD_012(W, 2, 0, 1)	if (((r_free >> 6) & 7) != S){\
r_free &= 077 | (S << 6);	UPD_AND_CL2(W, 2)}\
FD[W][0].bf.u32[3] = r_free;}

ZHGXN zhgxn;
ZHOU2 zhou2[10];
void ZHGXN::SetupFsol(int * grid0) {
	g0 = grid0;
	memset(fsol, 0, sizeof fsol);
	for (int i = 0; i < 81; i++) 	fsol[grid0[i]].Set_c(i);
	//for (int i = 0; i < 9; i++) fsol[i].Print3("sol");
}

int ZHOU2::GoZ2(int  fl) {
	zhgxn.nua = 0;
	if (_popcnt32(fl) != 2) {
		cout << "bug fl not 2 digits" << endl;
		return 1;// not valid fl
	}
	//memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	//memcpy(FD, zhoustart, sizeof FD);
	int n = 0;
	BF128 isfl;	isfl.SetAll_0();
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zhgxn.fsol[i];
			zhgxn.fsolw[n++] = zhgxn.fsol[i];
		}
	}
	cells_unsolved = isfl;
	FD[0][1] = isfl; FD[1][1] = isfl;
	isfl.bf.u32[3] = 0777777;
	FD[0][0] = isfl; FD[1][0] = isfl;
	//ImageCandidats();

	// start with guess on first cell
	uint32_t  cell = 10;
	{
		register uint32_t row1 = cells_unsolved.bf.u32[0], bit = 1, i = 0;
		for (; i < 9; i++, bit <<= 1)
			if (row1&bit) { cell = i; break; }
	}
	if (cell > 8)return 1; //this would be bug
	int xcell = C_To128[cell];
	zhou2[1] = *this;
	zhou2[1].Assign(0, cell, xcell);
	zhou2[1].ComputeNext();
	zhou2[1] = *this;
	zhou2[1].Assign(1, cell, xcell);
	zhou2[1].ComputeNext();
	return 0;
}
void ZHOU2::Assign(int digit, int cell, int xcell) {
	FD[digit][0] &= AssignMask_Digit[cell];
	cells_unsolved.Clear(xcell);
	if (digit)  FD[0][0].Clear(xcell);
	else FD[1][0].Clear(xcell);
}

int ZHOU2::ApplySingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	BF128 R1, R2;
	R1 = FD[0][0] | FD[1][0];
	R2 = FD[0][0] & FD[1][0];
	R1 -= R2;
	R1 &= cells_unsolved;
	if (R1.isEmpty()) {// no single store pairs
		zhgxn.cell_to_guess = R2.getFirstCell();
		return 0;
	}
	zh_g.single_applied = 1;
	int cell, xcell;
	while ((cell = R1.getFirstCell()) >= 0) {
		R1.Clear_c(cell);
		xcell = C_To128[cell];
		if (FD[0][0].On_c(cell))		Assign(0, cell, xcell);
		else if (FD[1][0].On_c(cell))	Assign(1, cell, xcell);
		else return 1;// conflict with previous assign
	}
	return 0;
}

int ZHOU2::Update() {
	register uint32_t Shrink = 1, r_free, B, A, S, last_assigned = 0, cur_assigned;
loop_upd:
	cur_assigned = last_assigned; last_assigned = 0;
	if (FD[1][0].bf.u32[3]) { UPD_ONE_DIGIT2(1) }
	if (FD[0][0].bf.u32[3]) { UPD_ONE_DIGIT2(0) }
exit_digits:
	if (last_assigned) goto loop_upd;// nothing to do in the next cycle
	return 1;
}
int ZHOU2::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU2::ComputeNext() {
	int ir = FullUpdate();
	if (!ir) return;// locked
	if (ir == 2) {//solved
		FD[0][0].bf.u32[3] = 0;
		if (FD[0][0] != zhgxn.fsolw[0]) {
			BF128 w0 = FD[0][0] - zhgxn.fsolw[0],
				w1 = FD[1][0] - zhgxn.fsolw[1];// this is the uaw0.Print3("ua");			
			zhgxn.tua[zhgxn.nua++] = w0 | w1;
			//(w0 | w1).Print3("added");
		}
		return;
	}
	Guess();// continue the process
}
void ZHOU2::Guess() {
	int cell = zhgxn.cell_to_guess,
		xcell = C_To128[cell];
	if (FD[0][0].On(xcell)) {
		ZHOU2 * mynext = (this + 1);
		*mynext = *this;
		mynext->Assign(0, cell, xcell);
		mynext->ComputeNext();
	}
	if (FD[1][0].On(xcell)) {
		ZHOU2 * mynext = (this + 1);
		*mynext = *this;
		mynext->Assign(1, cell, xcell);
		mynext->ComputeNext();
	}
}
void ZHOU2::ImageCandidats() {
	BF128  R2 = FD[0][0] & FD[1][0], R1 = FD[0][0] | FD[1][0];
	for (int i = 0; i < 9; i++) { // rows
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < 30; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			if (R1.Off(xcell)) cout << "-  ";
			else if (R2.On(xcell)) cout << "12 ";
			else if (FD[0][0].On(xcell))cout << "1  ";
			else cout << "2  ";
		} // end for j
		cout << endl;
	}
	cout << endl;

}

//====================== code for UA gen3 digits

#define UPD_AND_CL3(W,X)	last_assigned = W+1;cur_assigned = 0;\
B = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[X] &= B;\
FD[0][0].bf.u32[X] &= B; FD[1][0].bf.u32[X] &= B; FD[2][0].bf.u32[X] &= B;  }\
FD[W][0].bf.u32[X] = FD[W][1].bf.u32[X] = A;


#define UPD_ONE_DIGIT3(W) if (cur_assigned > W+1)goto exit_digits;\
	if (FD[W][0].bf.u64[0] != FD[W][1].bf.u64[0]\
|| FD[W][0].bf.u32[2] != FD[W][1].bf.u32[2]){\
r_free = FD[W][0].bf.u32[3];\
UPD_012(W, 0, 1, 2)	if ((r_free & 7) != S){\
r_free &= 0770 | S;	UPD_AND_CL3(W, 0)}\
UPD_012(W, 1, 0, 2)	if (((r_free >> 3) & 7) != S){\
r_free &= 0707 | (S << 3);	UPD_AND_CL3(W, 1)}\
UPD_012(W, 2, 0, 1)	if (((r_free >> 6) & 7) != S){\
r_free &= 077 | (S << 6);	UPD_AND_CL3(W, 2)}\
FD[W][0].bf.u32[3] = r_free;}

ZHOU3 zhou3[10];

int ZHOU3::GoZ3(int  fl) {
	if (_popcnt32(fl) != 3) return 1;// not valid fl
	//memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	//memcpy(FD, zhoustart, sizeof FD);
	int n = 0;
	BF128 isfl;	isfl.SetAll_0();
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zhgxn.fsol[i];
			zhgxn.digit_map[i] = n;
			zhgxn.fsolw[n++] = zhgxn.fsol[i];
		}
	}
	cells_unsolved = isfl;
	FD[0][1] = isfl; FD[1][1] = isfl; FD[2][1] = isfl;
	isfl.bf.u32[3] = 0777777;
	FD[0][0] = isfl; FD[1][0] = isfl; FD[2][0] = isfl;
	return 0;
	ImageCandidats();
	// start with guess on one cell box1  one cell box 4
	int  cell1, cell2;
	BF128 w = cells_unsolved;	w &= units3xBM[18];// box1

	cell1 = w.getFirstCell();
	w = cells_unsolved; w &= units3xBM[22];/// box4
	cell2 = w.getFirstCell();
	if (cell1 < 0 || cell2 < 0) return 1; // should be bug

	int xcell1 = C_To128[cell1], xcell2 = C_To128[cell2];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			zhou3[1] = *this;
			zhou3[1].Assign(i, cell1, xcell1);
			zhou3[1].Assign(j, cell2, xcell2);
			//zhou3[1].ImageCandidats();
			zhou3[1].ComputeNext();
		}
	return 0;
}

int ZHOU3::DoZ3(int * t, int nt) {// called in zhou3[1]
	zhgxn.nua = 0;
	*this = zhou3[0];
	for (int i = 0; i < nt; i++) {
		int cell = t[i], xcell = C_To128[cell],
			digit = zhgxn.g0[cell], dmap = zhgxn.digit_map[digit];
		Assign(dmap, cell, xcell);
	}
	//ImageCandidats();
	ComputeNext();
	return 0;
}

void ZHOU3::Assign(int digit, int cell, int xcell) {
	FD[digit][0] &= AssignMask_Digit[cell];
	cells_unsolved.Clear(xcell);
	if (!digit) { FD[1][0].Clear(xcell); FD[2][0].Clear(xcell); }
	else if (digit == 1) { FD[0][0].Clear(xcell); FD[2][0].Clear(xcell); }
	else { FD[0][0].Clear(xcell); FD[1][0].Clear(xcell); }
}

int ZHOU3::ApplySingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	BF128 R1, R2, R3;
	R1 = FD[0][0] | FD[1][0];
	R2 = FD[0][0] & FD[1][0];
	R3 = R2 & FD[2][0];	R2 |= (R1&FD[2][0]);	R1 |= FD[2][0];
	R1 -= R2;
	R2 -= R3;
	R1 &= cells_unsolved;
	if (R1.isEmpty()) {// no single store pairs
		if (R2.isNotEmpty())		zhgxn.cell_to_guess = R2.getFirstCell();
		else zhgxn.cell_to_guess = R3.getFirstCell();
		return 0;
	}
	zh_g.single_applied = 1;
	int cell, xcell;
	while ((cell = R1.getFirstCell()) >= 0) {
		R1.Clear_c(cell);
		xcell = C_To128[cell];
		if (FD[0][0].On_c(cell))		Assign(0, cell, xcell);
		else if (FD[1][0].On_c(cell))	Assign(1, cell, xcell);
		else if (FD[2][0].On_c(cell))	Assign(2, cell, xcell);
		else return 1;// conflict with previous assign
	}
	return 0;
}

int ZHOU3::Update() {
	register uint32_t Shrink = 1, r_free, B, A, S, last_assigned = 0, cur_assigned;
loop_upd:
	cur_assigned = last_assigned; last_assigned = 0;
	if (FD[2][0].bf.u32[3]) { UPD_ONE_DIGIT3(2) }
	if (FD[1][0].bf.u32[3]) { UPD_ONE_DIGIT3(1) }
	if (FD[0][0].bf.u32[3]) { UPD_ONE_DIGIT3(0) }
exit_digits:
	if (last_assigned) goto loop_upd;// nothing to do in the next cycle
	return 1;
}
int ZHOU3::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU3::ComputeNext() {
	int ir = FullUpdate();
	//cout << "back full ir=" << ir << endl;
	//ImageCandidats();
	if (!ir) return;// locked
	if (ir == 2) {//solved
		BF128 w0 = FD[0][0] - zhgxn.fsolw[0],
			w1 = FD[1][0] - zhgxn.fsolw[1],
			w2 = FD[2][0] - zhgxn.fsolw[2];// this is the ua;	
		//(w0 | w1 | w2).Print3("ua");
		if (w0.isEmpty() || w1.isEmpty() || w2.isEmpty())return;// already seen
		zhgxn.tua[zhgxn.nua++] = w0 | w1 | w2;
		return;
	}
	Guess();// continue the process
}
void ZHOU3::Guess() {
	int cell = zhgxn.cell_to_guess,
		xcell = C_To128[cell];
	if (FD[0][0].On(xcell)) {
		ZHOU3 * mynext = (this + 1);
		*mynext = *this;
		mynext->Assign(0, cell, xcell);
		mynext->ComputeNext();
	}
	if (FD[1][0].On(xcell)) {
		ZHOU3 * mynext = (this + 1);
		*mynext = *this;
		mynext->Assign(1, cell, xcell);
		mynext->ComputeNext();
	}
	if (FD[2][0].On(xcell)) {
		ZHOU3 * mynext = (this + 1);
		*mynext = *this;
		mynext->Assign(2, cell, xcell);
		mynext->ComputeNext();
	}
}
void ZHOU3::ImageCandidats() {
	BF128  R2 = FD[0][0] & FD[1][0], R1 = FD[0][0] | FD[1][0];
	BF128 R3 = R2 & FD[2][0];
	R2 |= R1 & FD[2][0];
	R1 |= FD[2][0];
	for (int i = 0; i < 9; i++) { // rows
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < 38; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			if (R1.Off(xcell)) cout << "-   ";
			else if (R3.On(xcell)) cout << "123 ";
			else if (R2.On(xcell)) {
				if (FD[0][0].Off(xcell))cout << "23  ";
				else if (FD[1][0].Off(xcell))cout << "13  ";
				else cout << "12  ";
			}
			else {
				if (FD[0][0].On(xcell))cout << "1   ";
				else if (FD[1][0].On(xcell))cout << "2   ";
				else cout << "3   ";
			}
		} // end for j
		cout << endl;
	}
	cout << endl;
}