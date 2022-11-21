/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
#include "main.h"   
#include "Zh1b2b.h"
extern uint64_t  p_cpt2g[70];

// row unsolved 6 bits per digit lookup table
uint64_t zh2b_t_runsolved[9] = { 077 , 077 << 6 , 077 << 12 , 077 << 18 , 077 << 24 ,
(uint64_t)077 << 32 ,(uint64_t)077 << 38 ,(uint64_t)077 << 44 ,(uint64_t)077 << 50 };
uint32_t zh2b_t_runsolvedshift[9] = { 0,6,12,18,24,32,38,44,50 };
//========================= global variable and working areas 
extern ZH_GLOBAL zh_g;
ZH2B zh2b_i, zh2b_i1,zh2b[40] ;
ZH2B_GLOBAL   zh2b_g;   // 2 bands 9 digits
ZH2GXN zh2gxn;
ZH2_2  zh2_2[5];
ZH2_3  zh2_3[10];
ZH2_4  zh2_4[20];
ZH2_5  zh2_5[30];

//============================= ZH_GLOBAL code
ZH2B_GLOBAL::ZH2B_GLOBAL(){
	zsol = 0; // no solution unless required buy the user
	nctlg =  0;

}

void ZH2B_GLOBAL::Init_g0(int* g0) {
	ndigits = 9;
	memcpy(puz0, g0, sizeof puz0);
	memset(fd_sols, 0, sizeof fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {// i column
		for (int j = 0; j < 3; j++) {//j row band1 or row band 2
			int cell = 9 * j + i, dig = puz0[cell];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = puz0[cell + 27];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
}
void ZH2B_GLOBAL::InitGangster(int * g0, int * g1) {
	uint64_t col64 = (uint64_t)Zhoucol | ((uint64_t)Zhoucol << 32);
	memcpy(fd_revised, fd_sols[1], sizeof fd_revised);
	for (int i = 0; i < 9; i++, col64 <<= 1) {
		if (g0[i] == g1[i])continue;
		int changes = g0[i] ^ g1[i]; // one added one cleared
		for (int d = 0, bit = 1; d < 9; d++, bit <<= 1) {// check digits
			if (!(changes & bit)) continue;
			if (g0[i] & bit)fd_revised[d] &= ~col64;
			else fd_revised[d] |= col64;
		}		
	}
}
uint64_t ZH2B_GLOBAL::BuildUaret(BF64 * wsol) {
	ua_ret = 0;
	for (int i = 0; i < 9; i++) {
		BF64 w = wsol[i] - fd_sols[0][i];
		ua_ret |= w.bf.u64;
	}
	return ua_ret;
}


//============= zh2B code for uas 2 bands and validity 2 bands puzzle

void ZH2B::Init_gang() {//init after zh1b_g InitGangster
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_revised, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}

#define UPDN(I,J)A=FD[I].bf.u32[J];\
Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
FD[I].bf.u32[1 - J] &= TblMaskSingle[S]; \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
CompFD[I].bf.u32[J] = FD[I].bf.u32[J] = A;

#define UPWCL(I,P,Q,R,T,U,V,W,X)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;

#define UPWCL2(I,P)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; wcl[P]&= cl;

#define UPWCL3(I,P,Q)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;

#define UPWCL4(I,P,Q,R)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;

#define UPWCL5(I,P,Q,R,T)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;


int ZH2B::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map, R4;// digits 12
	R2 |= R1 & Map; R1 |= Map;

	Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[4];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[5];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[6];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[7];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[8];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		//if (zh2b_g.diag) cout << Char2Xout(R1) << " apply R1 " << endl;
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 9; idig++) {
				if (map[idig] & bit) {// this is the digit
					int cell = From_128_To_81[res];
					Assign(idig, cell, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		R3 &= ~R4;
		uint32_t res, cell;
		if (zh2b_g.diag) 	cout<<Char2Xout(R2) << "R2  count "<< _popcnt64(R2) << endl;
		if (!R2) {
			if (!R3) R3 = R4;
			bitscanforward64(res, R3);
			zh2b_g.guess_xcell = res;
			return 0;
		}
		// try to get 2 cells or more
		while (bitscanforward64(res, R2)) {
			zh2b_g.guess_xcell = res;
			cell = From_128_To_81[res];
			uint64_t mask=R2 & cell_z3x[cell].u64[0];
			if (_popcnt64(mask)) return 0;
			if (_popcnt64(R2) < 3)return 0;
			R2 ^=( uint64_t)1<< res;
		}
		return 0;
	}
}
int ZH2B::FullUpdate() {
	//if (zh2b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}


//================================ZH2B code

inline int ZH2B::Seta(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	if (digit > 4) ddig += 2;// second bloc of 32 bits
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[8];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
inline void ZH2B::Assign(int digit, int cell, int xcell) {
	register uint64_t bit = (uint64_t)1 << xcell;
	FD[digit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.bf.u64 &= ~bit;
	int ddig = 6 * digit;
	if (digit > 4) ddig += 2;// second bloc of 32 bits
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	//if (digit != zh2b_g.puz0[cell])	cells_solved_false|= bit;
}

int ZH2B::Update(){
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink ){
	  Shrink = 0;
	  if (!rows_unsolved.bf.u32[0])goto digit5;

	  {register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
	  if (!(AR & 077))goto digit1;

//=digit 0
	  if (FD[0].bf.u32[0] ==  CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0,0)	if ((AR & 7) != S){
				AR &= 07777777770 | S;	UPWCL(0, 2,4,6,8,10,12,14,16)	}

digit0b:if (FD[0].bf.u32[1] ==CompFD[0].bf.u32[1])goto digit1;
		UPDN(0,1)	if (((AR >> 3) & 7) != S){
				AR &= 07777777707 | (S << 3);	UPWCL(1, 3,5,7,9,11,13,15,17)	}

digit1:	if (!(AR  & 07700))goto digit2;

		if (FD[1].bf.u32[0] ==	CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1,0)	if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);UPWCL(0, 0, 4, 6, 8, 10, 12, 14, 16)	}

digit1b:if (FD[1].bf.u32[1] ==	CompFD[1].bf.u32[1])goto digit2;
		UPDN(1,1)		if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9);UPWCL(1, 1, 5, 7, 9, 11, 13, 15, 17)	}

digit2:	if (!(AR  & 0770000))goto digit3;

		if (FD[2].bf.u32[0] ==	CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2,0)	if (((AR >> 12) & 7) != S){
			AR &= 07777707777 | (S << 12);	UPWCL(0, 0, 2, 6, 8, 10, 12, 14, 16)}

digit2b:if (FD[2].bf.u32[1] ==	CompFD[2].bf.u32[1])goto digit3;
		UPDN(2,1)	if (((AR >> 15) & 7) != S){
			AR &= 07777077777 | (S << 15);	UPWCL(1, 1, 3, 7, 9, 11, 13, 15, 17)	}

digit3: if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		  UPDN(3,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18); UPWCL(0, 0, 2, 4, 8, 10, 12, 14, 16)	 }

digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		  UPDN(3, 1)if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21); UPWCL(1, 1, 3, 5, 9, 11, 13, 15, 17)	}

digit4:if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] ==	CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S){
			AR &= 07077777777 | (S << 24);  UPWCL(0, 0, 2, 4, 6, 10, 12, 14, 16)	}

digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S){
			AR &= 0777777777 | (S << 27);  UPWCL(1, 1, 3, 5, 7, 11, 13, 15, 17)	}

end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


digit5:
	  if (!rows_unsolved.bf.u32[1])continue; // second lot  4 digits

	  {register unsigned int  AR = rows_unsolved.bf.u32[1];// valid for digits 5,6,7,8
	  if (!(AR & 077))goto digit6;

	  if (FD[5].bf.u32[0] ==	CompFD[5].bf.u32[0])goto digit5b;
		UPDN(5, 0) if ((AR & 7) != S){
			AR &= 07777777770 | S;	UPWCL(0, 0, 2, 4, 6, 8, 12, 14, 16)		}

digit5b:if (FD[5].bf.u32[1] == CompFD[5].bf.u32[1])goto digit6;
		UPDN(5, 1) 	if (((AR >> 3) & 7) != S){
			AR &= 07777777707 | (S << 3); UPWCL(1, 1, 3, 5, 7, 9, 13, 15, 17)	}

digit6:  if (!(AR & 07700))goto digit7;

		if (FD[6].bf.u32[0] ==  CompFD[6].bf.u32[0])goto digit6b;
		UPDN(6, 0) if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);	UPWCL(0, 0, 2, 4, 6, 8, 10, 14, 16)	  }

digit6b: if (FD[6].bf.u32[1] == CompFD[6].bf.u32[1])goto digit7;
		UPDN(6, 1) if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9); UPWCL(1, 1, 3, 5, 7, 9, 11, 15, 17)  }

digit7:   if (!(AR  & 0770000))goto digit8;

		if (FD[7].bf.u32[0] ==			  CompFD[7].bf.u32[0])goto digit7b;
		  UPDN(7, 0)if (((AR >> 12) & 7) != S){
			  AR &= 07777707777 | (S << 12);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 16)  }

digit7b:  if (FD[7].bf.u32[1] ==  CompFD[7].bf.u32[1])goto digit8;
		  UPDN(7, 1) if (((AR >> 15) & 7) != S){
			  AR &= 07777077777 | (S << 15);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 17)  }

digit8:  if (!(AR  & 077000000))goto end5678;

		  if (FD[8].bf.u32[0] ==  CompFD[8].bf.u32[0])goto digit8b;
		  UPDN(8,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 14)	  }

digit8b: if (FD[8].bf.u32[1] == CompFD[8].bf.u32[1])goto end5678;
		  UPDN(8,1)	  if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 15)  }

end5678:rows_unsolved.bf.u32[1] = AR;
	  }// end of validity for AR
  }// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
  return 1;
}
char * ZH2B::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	int tdig[9];// build the table of digits
	register int A = rows_unsolved.bf.u32[0];
	tdig[0] = A & 077; A >>= 6; tdig[1] = A & 077; A >>= 6; tdig[2] = A & 077; A >>= 6;
	tdig[3] = A & 077; A >>= 6; tdig[4] = A & 077;
	A = rows_unsolved.bf.u32[1];
	tdig[5] = A & 077; A >>= 6; tdig[6] = A & 077; A >>= 6;
	tdig[7] = A & 077; A >>= 6; tdig[8] = A & 077;

	for (int idig = 0; idig < 9; idig++) {// one digit
		int arowsj = tdig[idig];// 6 rows
		if (arowsj == 077) continue;
		for (int ib = 0; ib < 2; ib++) {// 3 blocs per digit
			int arows = (arowsj >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = FD[idig].bf.u32[ib];
			for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				uint32_t  irow;
				bitscanforward(irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = idig + '1';
			}
		}
	}
	return zs;
}
void ZH2B::Init_std_bands() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
inline void ZH2B::InitTclues(uint32_t * tclues, int n) {
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int icell = 0; icell < n; icell++) {
		int cell = tclues[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map

		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
}


void ZH2B::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh2b_g.nsol << " unsolved=" << Unsolved_Count() << endl;
	//	cout << zh1b_g.out27 << " band1 "<<endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_6 digits " << endl;
	if (!all) return;

	cout << "map per digit bands 2 3" << endl;
	for (int ib = 0; ib < 2; ib++) {
		for (int ir = 0; ir < 3; ir++) {
			for (int idig = 0; idig < 9; idig++) {
				unsigned vf = FD[idig].bf.u32[ib];
				unsigned wr = (vf >> (9 * ir)) & 0x1ff;
				for (int k = 0; k < 9; k++) {
					if (wr & (1 << k))		cout << idig + 1;
					else 		cout << ".";
					if (k == 2 || k == 5) 	cout << " ";
				}
				cout << "  ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
int ZH2B::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (int i = 0; i < 9; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " 
		<< rows_unsolved.Count() << " non assigned cells"  << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres imposés
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}


void ZH2B::InitBands12(int * g0) {
	zh2b_g.ndigits = 9;
	memcpy(zh2b_g.puz0, g0, sizeof zh2b_g.puz0);
	memset(zh2b_g.fd_sols, 0, sizeof zh2b_g.fd_sols);
	// build sol per digit and pm per digit at start
	for (int i = 0; i < 9; i++) {// i column
		for (int j = 0; j < 3; j++) {//j row band1 or row band 2
			int cell = 9 * j + i, dig = zh2b_g.puz0[cell];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = zh2b_g.puz0[cell + 27];
			zh2b_g.fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			zh2b_g.fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			zh2b_g.fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	//ImageCandidats();

}

uint64_t ZH2B::IsValid(uint32_t * tclues, int n,int onlyone) {
	*this = zh2b[0];
	InitTclues(tclues, n);
	zh2gxn.nua = 0;	zh2gxn.uamin = 100;
	zh2gxn.onlyone = onlyone;
	zh2b_g.go_back = 0;
	ComputeNext();
	return zh2gxn.nua;
}
uint64_t ZH2B::IsValid(uint64_t bf54,  int onlyone) {
	*this = zh2b[0];
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	register uint64_t B = bf54;
	int cell;
	while (bitscanforward64(cell, B)) {
		B ^= (uint64_t)1 << cell;
		int  digit = zh2b_g.puz0[cell], xcell = C_To128[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
	zh2gxn.nua = 0;	zh2gxn.uamin = 100;
	zh2gxn.onlyone = onlyone;
	zh2b_g.go_back = 0;
	ComputeNext();
	return zh2gxn.nua;
}

inline void ZH2B::ComputeNext() {
	if (zh2b_g.go_back) return;
	int ir = FullUpdate();
	//cout << "back full update ir="<<ir	<< endl;		
	//ImageCandidats();
	if (ir == 1)GuessValidB12();
	else if (ir == 2) {// solved 
		zh2gxn.uaw = 0;
		for (int i = 0; i < 9; i++) {
			BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
			zh2gxn.uaw |= w.bf.u64;
		}
		if (zh2gxn.uaw) {
			uint64_t cc = _popcnt64(zh2gxn.uaw);
			if (cc < zh2gxn.uamin)zh2gxn.uamin = cc;
			if (cc > zh2gxn.uamin) return;
			zh2gxn.tua[zh2gxn.nua++] = zh2gxn.uaw;
			//cout << Char2Xout(zh2gxn.uaw) << " uaw "
				//<< cc << " " << zh2gxn.nua << endl;
			if (zh2gxn.onlyone)zh2b_g.go_back = 1;
		}
	}
}
void ZH2B::GuessValidB12() {// 
	if (zh2b_g.go_back) return;
	uint32_t xcell = zh2b_g.guess_xcell, cell, digit;
	cell = From_128_To_81[xcell];
	digit = zh2b_g.puz0[cell];
	uint64_t bit = (uint64_t)1 << xcell;
	// true first if possible
	if (FD[digit].bf.u64 & bit) {
		//cout << digit+1 << cellsFixedData[cell].pt << " try true" << endl;
		//cout << "Guess okr " << digit + 1 << cellsFixedData[cell].pt << endl;
		ZH2B * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
		if (zh2b_g.go_back) return;
	}
	// then false 
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig].bf.u64 & bit) {
			//cout << "Guess nokr " << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2B * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
			if (zh2b_g.go_back) return;
		}
	}

}

//______________ uas collector using uas

void ZH2B::InitB1245() {
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 6; col < 9; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	Do4bGo();

}
void ZH2B::InitB1346() {
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 3; col < 6; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	Do4bGo();
}
void ZH2B::InitB2356() {
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign box 3 6
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	for (int row = 0; row < 6; row++) for (int col = 0; col < 3; col++) {
		int cell = 9 * row + col, xcell = C_To128[cell],
			digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	Do4bGo();
}
void ZH2B::InitBf(uint64_t bf) {
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign bf
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	int xcell;
	while (bitscanforward64(xcell, bf)) {
		bf ^= (uint64_t)1 << xcell;
		int cell = From_128_To_81[xcell], digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
}
void ZH2B::InitBfG2(uint64_t bf, int c1, int d1, int c2, int d2) {
	memcpy(this, zh2b_start, sizeof zh2b_start);
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// assign bf
	memset(zh2b_g.Digit_cell_Assigned_step, 0, sizeof zh2b_g.Digit_cell_Assigned_step);
	int xcell;
	register uint64_t Bf = bf;
	while (bitscanforward64(xcell, Bf)) {// must be 3 digits all cells 
		Bf ^= (uint64_t)1 << xcell;
		int cell = From_128_To_81[xcell], digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_step[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_step[i];
	int bit12 = (1 << d1) | (1 << d2);
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i],col= C_col[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		if (col == c1|| col == c2) {
			FD[d1].bf.u64 ^= bit; 
			FD[d2].bf.u64 ^= bit;
		}
	}
	//cout << Char2Xout(bf)
	//	<< " cols " << c1 + 1 << c2 + 1 << " digs " << d1 + 1 << d2 + 1 << endl;
	//ImageCandidats();
}

int ZH2B::Do4bGo() {
	zh2b_g.limadd = 25;
	zh2b_g.test4b = zh2b_g.ntest4b = 0;
	zh2gxn.nua = 0;
	zh2b_g.go_back = 0;
	ComputeNextb();
	return 0;
}
int ZH2B::Dob(int lim) {
	zh2b_g.limadd = lim;
	zh2b_g.test4b = 1;
	zh2gxn.nua = zh2b_g.ntest4b = 0;
	zh2b_g.go_back = 0;
	ComputeNextb();
	return zh2b_g.ntest4b;
}

int ZH2B::FullUpdateb() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCellsb())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2B::ApplySingleOrEmptyCellsb() {
	zh2b_g.single_applied = 0;
	uint64_t* map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map, R4, R5, R6;// digits 12
	R2 |= R1 & Map; R1 |= Map;

	Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	Map = map[4];  R5 = R4 & Map;
	R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[5]; R6 = R5 & Map; R5 |= R4 & Map;
	R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[6]; R6 |= R5 & Map; R5 |= R4 & Map;
	R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[7]; R6 |= R5 & Map; R5 |= R4 & Map;
	R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[8]; R6 |= R5 & Map; R5 |= R4 & Map;
	R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		//if (zh2b_g.diag) cout << Char2Xout(R1) << " apply R1 " << endl;
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 9; idig++) {
				if (map[idig] & bit) {// this is the digit
					int cell = From_128_To_81[res];
					Assign(idig, cell, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;	R3 &= ~R4; R4 &= ~R5; R5 &= ~R6;
		if (!R2) { R2 = R3; {if (!R2) R2 = R4; {if (!R2) R2 = R5; }} }
		uint64_t sok = 0; // zh2gxn.unsolved_field;
		for (int i = 0; i < 9; i++) {
			sok |= (FD[i].bf.u64 & zh2gxn.fsol[i]);
		}
		uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
		uint64_t nok = 54 - _popcnt64(cells_unsolved.bf.u64) - _popcnt64(sok2);
		if (nok >= zh2b_g.limadd) return -1;// already too big
		sok &= cells_unsolved.bf.u64;
		int nu = *zh2gxn.nknownuas;
		for (int iu = 0; iu < nu; iu++) {
			register uint64_t U = zh2gxn.knownuas[iu];
			if (U & sok2) continue;
			U &= sok;
			if (!U)return -1; // dead branch
			if (U & R2)U &= R2;
			else if (U & R3)U &= R3;
			else if (U & R4)U &= R4;
			else if (U & R5)U &= R5;
			else U &= R6;
			//if(zh2b_g.test4b == 1)
			//	cout << Char2Xout(U) << "U to fill iu=" <<iu<< endl;
			bitscanforward64(zh2b_g.guess_xcell, U);
			return 0;
		}
		bitscanforward64(zh2b_g.guess_xcell, R2);
		return 0;
	}
}
inline void ZH2B::ComputeNextb() {
	if (zh2b_g.go_back) return;
	int ir = FullUpdateb();
	//cout << "back full ir=" << ir << endl;
	if (!ir) return;
	//cells_solved_false = 0;
	for (int i = 0; i < 9; i++) {
		BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
		//cells_solved_false |= w.bf.u64;
	}
	//cells_solved_false &= ~cells_unsolved.bf.u64;
	if (ir == 1)Guessb();
	else if (ir == 2) {// solved 
		zh2b_g.ntest4b++;
		uint64_t ww = 0;
		for (int i = 0; i < 9; i++) {
			BF64 w = FD[i] - zh2b_g.fd_sols[0][i];
			ww |= w.bf.u64;
		}
		if (ww) {
			uint64_t cc = _popcnt64(ww);
			//cout << Char2Xout(ww) << " seen  cc" << cc << endl;
			if (cc < 12) return;
			if (cc >= zh2b_g.limadd) return;
			// check no ua false
			int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 9; i++) {
					sok |= (FD[i].bf.u64 & zh2gxn.fsol[i]);
				}
				//cout << Char2Xout(sok) << " sok" << endl;
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U & sok)) return;;
				}
			}

			zh2gxn.tua[zh2gxn.nua++] = ww;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++] = ww;
			//cout << Char2Xout(ww) << "4b added zh2gxn.nua="<< zh2gxn.nua 
				//<< " cc="<<cc<< endl;
			if (zh2gxn.nua > 80)zh2b_g.go_back = 1;
		}
	}
}
void ZH2B::Guessb() {// next free ua 
	if (zh2b_g.go_back) return;
	//cout << Char2Xout(cells_solved_false) << "guessb false index"
	//	<<this-&zh2b[0] << endl;
	//cout << Char2Xout(cells_unsolved.bf.u64)<<" unsolved" << endl;
	//ImageCandidats();
	uint32_t xcell = zh2b_g.guess_xcell, cell, digit;
	cell = From_128_To_81[xcell];
	digit = zh2b_g.puz0[cell];
	uint64_t bit = (uint64_t)1 << xcell;
	// true first if possible
	if (FD[digit].bf.u64 & bit) {
		//cout << "guess ok" <<digit+1<< cellsFixedData[cell].pt << endl;
		ZH2B* mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNextb();
	}
	// then false 
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (FD[idig].bf.u64 & bit) {
			//cout << "guess nok" << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2B* mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNextb();
		}
	}
}


//____________ uas collector 

void ZH2GXN::SetupFsol(int * grid0) {
	g0 = grid0;
	memset(fsol, 0, sizeof fsol);
	memset(gangsters, 0, sizeof gangsters);
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		fsol[g0[i]] |= bit;
		gangsters[C_col[i]] |= 1 << g0[i];
	}
}
//________________ 2 digits guas2 only 

void ZH2_2::GoZ2A(int fl) {
	int n = 0;
	zh2gxn.nua = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl & bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);

}

int ZH2_2::GoZ2D(int  fl) {

	*zh2gxn.nknownuas = 0;// be sure to start with no ua
	zh2gxn.nkguas = 0;
	GoZ2A(fl);// start shared with gangsters g2
	// init pm using gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		for (int idig = 0; idig < 2; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	diag = 0;
	//ImageCandidats();
	//return -1;
	if (!FullUpdate()) return -1;
	Guess();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_2::GoZ2G2(int fl, int c1, int d1, int c2, int d2) {
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl & bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
	// init pm using revised gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	int bit12 = (1 << d1) | (1 << d2);
	gx[c1] ^= bit12;// must do remove one add the other
	gx[c2] ^= bit12;

	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		for (int idig = 0; idig < 2; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	zh2gxn.nua = 0;

	diag = 0;
	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate ua size 18
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 2; i++) 
			w |= FD[i] - zh2gxn.fsolw[i];
		zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
		return 1;
	}
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
void ZH2_2::DoZ2Go() {
	zh2gxn.nua = 0;
	Guess();
}
inline void ZH2_2::Assign(int rdigit, int cell, int xcell) {
	FD[rdigit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.Clear(xcell);
	int ddig = 6 * rdigit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
}
int ZH2_2::Seta(int rdigit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell];
	if (FD[rdigit].Off(xcell)) return 1; // not valid
	Assign(rdigit, cell, xcell);
	if(rdigit)FD[0].Clear(xcell);
	else FD[1].Clear(xcell);
	return 0;
}
int ZH2_2::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2_2::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t* map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]);// digits 12
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 4; idig++) {
				if (map[idig] & bit) {// this is the digit
					Seta(idig, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		// setup ok cells
		uint64_t sok = 0; // zh2gxn.unsolved_field;
		for (int i = 0; i < 2; i++) {
			sok |= (FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64);
		}
		uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
		sok &= cells_unsolved.bf.u64;
		int nu = *zh2gxn.nknownuas;
		// first non hit ua to solve
		for (int iu = 0; iu < nu; iu++) {
			register uint64_t U = zh2gxn.knownuas[iu];
			if (diag)cout << Char2Xout(U) << " U iu=" << iu << endl;
			if (U & sok2) continue;
			U &= sok;
			if (!U)return -1; // dead branch
			bitscanforward64(zh2b_g.guess_xcell, U);
			return 0;
		}
		// all uas solved, first cell  
		bitscanforward64(zh2b_g.guess_xcell, R2);
		return 0;
	}
}
void ZH2_2::Guess() {
	uint32_t xcell = zh2b_g.guess_xcell, cell;
	cell = From_128_To_81[xcell];
	uint64_t bit = (uint64_t)1 << xcell;
	//ImageCandidats();
	int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
	if (FD[digit].bf.u64 & bit) {
		ZH2_2* mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
	}
	for (int idig = 0; idig < 3; idig++) {
		if (idig == digit) continue;
		if (FD[idig].bf.u64 & bit) {
			ZH2_2* mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
		}
	}
}
void ZH2_2::ComputeNext() {
	int ir = FullUpdate();
	if (ir == 1)Guess();
	else if (ir == 2) {// solved 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 2; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		uint64_t cc = _popcnt64(w.bf.u64);
		if (w.bf.u64) {
			// check no ua false
			int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 2; i++) {
					sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
				}
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U & sok)) return;;
				}
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
				= w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " added" << endl;
		}
	}
}
int ZH2_2::Update() {
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, * wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])break;

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL2(0, 2)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL2(1, 3)
		}

	digit1:	if (!(AR & 07700))goto end01;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL2(0, 0)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto end01;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL2(1, 1)
		}

	end01: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


	}// end while

	return 1;
}
void ZH2_2::ImageCandidats() {
	BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];

	for (int i = 0; i < 6; i++) { // rows
		if ((i == 3)) {
			for (int ix = 0; ix < 35; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			uint64_t	bit = (uint64_t)1 << xcell;
			if (!(R1.bf.u64 & bit)) cout << "-  ";
			else {
				for (int i = 0; i < 2; i++)
					if (FD[i].bf.u64 & bit)cout << i + 1;
				if (R2.bf.u64 & bit)cout << " ";
				else cout << "  ";
			}

		} // end for j
		cout << endl;
	}
}

//________________ 3 digits 
void ZH2_3::GoZ3A(int fl) {
	int n = 0;
	zh2gxn.nua = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
int ZH2_3::GoZ3(int  fl) {
	if (_popcnt32(fl) != 3) {
		cout << "bug fl not 3 digits" << endl;
		return -1;// not valid fl
	}
	*zh2gxn.nknownuas = 0;// be sure to start with no ua
	zh2gxn.nkguas = 0;
	GoZ3A(fl);// start shared with gangsters g2
	// init pm using gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 3; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	diag = 0;
	if(!FullUpdate()) return -1;
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_3::GoZ3G2(int fl, int c1, int d1, int c2, int d2) {
	if (_popcnt32(fl) != 3) {
		cout << "bug fl not 3 digits" << endl;
		return -1;// not valid fl
	}
	GoZ3A(fl);// start shared with gangsters g2

	// init pm using revised gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	int bit12 = (1 << d1)|(1 << d2);
	gx[c1] ^= bit12;// must do remove one add the other
	gx[c2] ^= bit12;

	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 3; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 0777777;//3*6 bits

	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate gua 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 3; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		register uint64_t uan = w.bf.u64;
		if (zh2gxn.nkguas) {// likely redundant
			for(uint32_t i=0;i< zh2gxn.nkguas;i++)
				if(!(~zh2gxn.knownuas[i] & (uan))) return 1;
		}
		zh2gxn.tua[zh2gxn.nua++] = uan;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_3::GoZ3G3(int fl, int* gx, int debug ) {
	if (_popcnt32(fl) != 3) {
		cout << "bug fl not 3 digits" << endl;
		return -1;// not valid fl
	}
	GoZ3A(fl);// start 
	// init pm using revised gangster
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		for (int idig = 0; idig < 3; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}

	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	//zh2gxn.nua = 0;
	if (debug) ImageCandidats();
	int ir = FullUpdate();
	if (debug) ImageCandidats();
	if (!ir) return -1;
	if (ir == 2) {// immediate gua
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 3; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		register uint64_t uan = w.bf.u64;
		if (zh2gxn.nkguas) {// likely redundant
			for (uint32_t i = 0; i < zh2gxn.nkguas; i++)
				if (!(~zh2gxn.knownuas[i] & (uan))) return 1;
		}
		zh2gxn.tua[zh2gxn.nua++] = uan;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_3::DoZ3Go( int debug) {
	if (debug) cout << "DoZ3Go with debug" << endl;
	zh2gxn.nua = 0;
	diag = debug;
	Guess();
	return 0;
}
inline void ZH2_3::Assign(int rdigit, int cell, int xcell) {
	FD[rdigit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.Clear(xcell);
	int ddig = 6 * rdigit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
}
int ZH2_3::Seta(int rdigit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell];
	if (FD[rdigit].Off(xcell)) return 1; // not valid
	Assign(rdigit, cell, xcell);
	BF64 *Fd = &FD[rdigit];
	BF64 * RF = &FD[2];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
int ZH2_3::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2_3::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map;// digits 12
	R2 |= R1 & Map; R1 |= Map;
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 4; idig++) {
				if (map[idig] & bit) {// this is the digit
					Seta(idig, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		if (!R2) R2 = R3;
		zh2gxn.rx = R2;
		// setup ok cells
		uint64_t sok = 0; // zh2gxn.unsolved_field;
		for (int i = 0; i < 3; i++) {
			sok |= (FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64);
		}
		uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
		sok &= cells_unsolved.bf.u64;
		int nu = *zh2gxn.nknownuas;
		if (diag) {
			cout << Char2Xout(sok) << " sok" << endl;
			cout << Char2Xout(sok2) << " sok2" << endl;
		}
		// first non hit ua to solve
		for (int iu = 0; iu < nu; iu++) {
			register uint64_t U = zh2gxn.knownuas[iu];
			if (diag)cout << Char2Xout(U) << " U iu="<<iu << endl;
			if (U & sok2) continue;
			U &= sok;
			if (!U)return -1; // dead branch
			if (U & R2) U &= R2;
			else U &= R3;
			bitscanforward64(zh2b_g.guess_xcell, U);
			return 0;
		}
		// all uas solved, first cell in rx
		bitscanforward64(zh2b_g.guess_xcell, R2);
		return 0;
	}
}

void ZH2_3::Guess() {
	uint32_t xcell = zh2b_g.guess_xcell, cell;
	cell = From_128_To_81[xcell];
	if (diag) {
		cout << "guess " << cell << " ";
		ImageCandidats();
	}
	uint64_t bit = (uint64_t)1 << xcell;
	int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
	if (FD[digit].bf.u64 & bit) {
		ZH2_3 * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
	}
	for (int idig = 0; idig < 3; idig++) {
		if (idig == digit) continue;
		if (FD[idig].bf.u64 & bit) {
			ZH2_3 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
		}
	}
}
void ZH2_3::ComputeNext() {
	if (diag) {
		cout << "compute next " ;
		ImageCandidats();
	}
	int ir = FullUpdate();
	if (ir == 1)Guess();
	else if (ir == 2) {// solved 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 3; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		uint64_t cc = _popcnt64(w.bf.u64);
		register uint64_t uan = w.bf.u64;
		int nu = *zh2gxn.nknownuas;
		if (diag) {
			cout << Char2Xout(uan) << " false nu" << nu << endl;
			ImageCandidats();
		}
		if (uan)	 {
			// check no ua false
			//int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 3; i++) {
					sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
				}
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu]
						& BIT_SET_2X;
					if (!(U &sok)) return;
					if(U == uan) return;
				}
			}
			zh2gxn.tua[zh2gxn.nua++] = uan;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
				= uan;
			if (diag)cout << Char2Xout(uan) << " added " << endl;
		} 
	}
}
int ZH2_3::Update() {
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])break;

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL3(0, 2, 4)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL3(1, 3, 5)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL3(0, 0, 4)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL3(1, 1, 5)
		}

	digit2:	if (!(AR & 0770000))goto end01234;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL3(0, 0, 2)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto end01234;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL3(1, 1, 3)
		}


	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


	}// end while

	return 1;
}
void ZH2_3::ImageCandidats() {
	BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
	BF64 R3 = R2 & FD[2];
	R2 |= R1 & FD[2];	R1 |= FD[2];

	for (int i = 0; i < 6; i++) { // rows
		if ((i == 3)) {
			for (int ix = 0; ix < 45; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			uint64_t	bit = (uint64_t)1 << xcell;
			if (!(R1.bf.u64&bit)) cout << "-   ";
			else if (R3.bf.u64&bit) cout << "123 ";
			else {
				for (int i = 0; i < 3; i++)
					if (FD[i].bf.u64&bit)cout << i + 1;
				if (R2.bf.u64&bit)cout << "  ";
				else cout << "   ";
			}

		} // end for j
		cout << endl;
	}
	//cout << endl;

}


//________________ 4 digits 
void ZH2_4::GoZ4A(int  fl) {
	zh2gxn.nua = 0;
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);

}
int ZH2_4::GoZ4(int  fl) {
	if (_popcnt32(fl) != 4) {
		cout << "bug fl not 4 digits" << endl;
		return -1;// not valid fl
	}
	GoZ4A( fl);
	// init pm using gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 4; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 077777777;//4*6 bits
	FullUpdate();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;;
}
int ZH2_4::GoZ4G2(int fl, int c1, int d1, int c2, int d2) {
	if (_popcnt32(fl) != 4) {
		cout << "bug fl not 4 digits" << endl;
		return -1;// not valid fl
	}
	GoZ4A(fl);// start shared with gangsters g2

	// init pm using revised gangster
	uint32_t gx[9];
	for (int i = 0; i < 9; i++) {
		gx[i] = zh2gxn.gangsters[i] & fl;
	}
	int bit12 = (1 << d1) | (1 << d2);
	gx[c1] ^= bit12;// must do remove one add the other
	gx[c2] ^= bit12;

	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64&bit))continue;
		for (int idig = 0; idig < 4; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}
	rows_unsolved.bf.u64 = 077777777;//4*6 bits


	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate ua size 18
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 4; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
		//cout << Char2Xout(w.bf.u64) << " seen " 
		//	<<_popcnt64(w.bf.u64 )<< endl;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_4::GoZ4G3(int fl, int* gx) {
	if (_popcnt32(fl) != 4) {
		cout << "bug fl not 4 digits" << endl;
		return -1;// not valid fl
	}
	GoZ4A(fl);// start 
	// init pm using revised gangster
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		for (int idig = 0; idig < 4; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}

	rows_unsolved.bf.u64 = 0777777;//3*6 bits

	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate gua
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 4; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		register uint64_t uan = w.bf.u64;
		if (zh2gxn.nkguas) {// likely redundant
			for (uint32_t i = 0; i < zh2gxn.nkguas; i++)
				if (!(~zh2gxn.knownuas[i] & (uan))) return 1;
		}
		zh2gxn.tua[zh2gxn.nua++] = uan;
		return 1;
	}
	//cout << "after update" << endl;
	//ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}
int ZH2_4::DoZ4Go() {// called in zhou3[1]
	zh2gxn.nua = 0;
	/*
	if (*zh2gxn.nknownuas) {
		int ir = GetNetUaCell();
		if (ir < 0)return -1; // dead branch
		zh2b_g.guess_xcell = ir; // dead branch
	}
	*/
	Guess();
	return 0;
}
inline void ZH2_4::Assign(int rdigit, int cell, int xcell) {
	FD[rdigit] &= AssignMask_Digit[cell].u64[0];
	cells_unsolved.Clear(xcell);
	int ddig = 6 * rdigit;
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
}
int ZH2_4::Seta(int rdigit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell];
	if (FD[rdigit].Off(xcell)) return 1; // not valid
	Assign(rdigit, cell, xcell);
	BF64 *Fd = &FD[rdigit];
	BF64 * RF = &FD[3];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}
int ZH2_4::FullUpdate() {
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())	return 0; // locked 
		if (zh2b_g.single_applied) 			continue;
		break;
	}
	return 1;
}
int ZH2_4::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	uint64_t * map = &FD[0].bf.u64;
	uint64_t unsolved = cells_unsolved.bf.u64;
	register uint64_t R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map = map[2],
		R3 = R2 & Map;// digits 12
	R2 |= R1 & Map; R1 |= Map;
	Map = map[3]; 	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // these are new singles	
	if (R1) {
		zh2b_g.single_applied = 1;
		while (R1) {// usually a very small number of cells to assign
			uint32_t res;
			if (!bitscanforward64(res, R1)) break;
			uint64_t bit = (uint64_t)1 << res; // switch to the bit value
			R1 &= ~bit;  // clear the bit
			for (int idig = 0; idig < 4; idig++) {
				if (map[idig] & bit) {// this is the digit
					Seta(idig, res);
					goto nextr1;// finished for that cell
				}
			}
			return 1; //conflict with a previous cell assugn
		nextr1: {}
		}
		return 0;
	}
	else {
		R2 &= ~R3;
		if (!R2) R2 = R3;
		// setup ok cells
		uint64_t sok = 0; // zh2gxn.unsolved_field;
		for (int i = 0; i < 4; i++) {
			sok |= (FD[i] & zh2gxn.fsolw[i]).bf.u64;
		}
		uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
		sok &= cells_unsolved.bf.u64;
		int nu = *zh2gxn.nknownuas;
		// first non hit ua to solve
		for (int iu = 0; iu < nu; iu++) {
			register uint64_t U = zh2gxn.knownuas[iu];
			if (U & sok2) continue;
			U &= sok;
			if (!U)return -1; // dead branch
			if (U & R2) U &= R2;
			else U &= R3;
			bitscanforward64(zh2b_g.guess_xcell, U);
			return 0;
		}
		// all uas solved, first cell in rx
		bitscanforward64(zh2b_g.guess_xcell, R2);
		return 0;

	}
}
void ZH2_4::Guess() {
	uint32_t xcell = zh2b_g.guess_xcell, cell;
	cell = From_128_To_81[xcell];
	uint64_t bit = (uint64_t)1 << xcell;
	//cout << "guess" << cellsFixedData[cell].pt << endl;
	int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
	if (FD[digit].bf.u64 & bit) {
		//cout << "guess ok" <<digit+1<< cellsFixedData[cell].pt << endl;
		ZH2_4 * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Seta(digit, xcell);
		mynext->ComputeNext();
	}
	for (int idig = 0; idig < 4; idig++) {
		if (FD[idig].bf.u64 & bit) {
			if (idig == digit) continue;
			//cout << "guess nok" << idig + 1 << cellsFixedData[cell].pt << endl;
			ZH2_4 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(idig, xcell);
			mynext->ComputeNext();
		}
	}
}
void ZH2_4::ComputeNext() {
	int ir = FullUpdate();
	//cout << "compnext ir=" << ir << endl;
	//ImageCandidats();
	if (ir == 1)Guess();
	else if (ir == 2) {// solved 
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 4; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		uint64_t cc = _popcnt64(w.bf.u64);
		if (w.bf.u64)	 {
			// check no ua false
			//cout << Char2Xout(w.bf.u64) << " seen" << endl;
			int nu = *zh2gxn.nknownuas;
			if (nu) {
				uint64_t sok = 0; // zh2gxn.unsolved_field;
				for (int i = 0; i < 4; i++) {
					sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
				}
				for (int iu = 0; iu < nu; iu++) {
					register uint64_t U = zh2gxn.knownuas[iu];
					if (!(U &sok)) return;;
				}
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
				= w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " added zh2gxn.nua="<< zh2gxn.nua << endl;
		}
	}
}
int ZH2_4::Update() {
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		if (!rows_unsolved.bf.u32[0])break;

		{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
		if (!(AR & 077))goto digit1;

		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL4(0, 2, 4, 6)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL4(1, 3, 5, 7)
		}

	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL4(0, 0, 4, 6)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL4(1, 1, 5, 7)
		}

	digit2:	if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL4(0, 0, 2, 6)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL4(1, 1, 3, 7)
		}

	digit3: if (!(AR & 077000000))goto end01234;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL4(0, 0, 2, 4)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto end01234;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL4(1, 1, 3, 5)
		}


	end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


	}// end while

	return 1;
}
void ZH2_4::ImageCandidats() {
	BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
	BF64 R3 = R2 & FD[2];
	R2 |= R1 & FD[2];	R1 |= FD[2];
	BF64 R4 = R3 & FD[3];
	R3 |= R2 & FD[3]; R2 |= R1 & FD[3];	R1 |= FD[3];

	for (int i = 0; i < 6; i++) { // rows
		if ((i == 3)) {
			for (int ix = 0; ix < 45; ix++)       cout << (char)'-';
			cout << endl;
		}
		for (int j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, xcell = C_To128[cell];
			uint64_t	bit = (uint64_t)1 << xcell;
			if (!(R1.bf.u64&bit)) cout << "-    ";
			else if (R4.bf.u64&bit) cout << "1234 ";
			else {
				for (int i = 0; i < 4; i++)
					if (FD[i].bf.u64&bit)cout << i + 1;
				if (R3.bf.u64&bit)cout << "  ";
				else if (R2.bf.u64&bit)cout << "   ";
				else cout << "    ";
			}

		} // end for j
		cout << endl;
	}
	cout << endl;

}

//_______________________________  5 digits
void ZH2_5::GoZ5A(int  fl) {
	zh2gxn.nua = 0;
	int n = 0;
	uint64_t isfl = 0;
	for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
		if (fl&bit) {
			isfl |= zh2gxn.fsol[i];
			zh2gxn.maptodigit[n] = i;
			zh2gxn.fsolw[n].bf.u64 = zh2gxn.fsol[i];
			zh2gxn.digit_map[i] = n++;
		}
	}
	cells_unsolved.bf.u64 = isfl;
	zh2gxn.unsolved_field = isfl;
	memset(FD, 0, sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
int ZH2_5::GoZ5(int  fl) {
		GoZ5A(fl);
		uint32_t gx[9];
		for (int i = 0; i < 9; i++) 
			gx[i] = zh2gxn.gangsters[i] & fl;		
		for (int i = 0; i < 54; i++) {
			int xi = C_To128[i];
			uint64_t bit = (uint64_t)1 << xi;
			if (!(cells_unsolved.bf.u64&bit))continue;
			for (int idig = 0; idig < 5; idig++) {
				uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
				if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
			}
		}
		rows_unsolved.bf.u64 = 07777777777;//5*6 bits
		uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
		return (int)cc;;
	}
int ZH2_5::GoZ5G2(int fl, int c1, int d1, int c2, int d2) {
		if (_popcnt32(fl) != 5) {
			cout << "bug fl not 5 digits" << endl;
			return -1;// not valid fl
		}
		GoZ5A(fl);
		uint32_t gx[9];
		for (int i = 0; i < 9; i++) {
			gx[i] = zh2gxn.gangsters[i] & fl;
		}
		int bit12 = (1 << d1) | (1 << d2);
		gx[c1] ^= bit12;// must do remove one add the other
		gx[c2] ^= bit12;
		for (int i = 0; i < 54; i++) {
			int xi = C_To128[i];
			uint64_t bit = (uint64_t)1 << xi;
			if (!(cells_unsolved.bf.u64&bit))continue;
			for (int idig = 0; idig < 5; idig++) {
				uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
				if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
			}
		}
		rows_unsolved.bf.u64 = 07777777777;//5*6 bits


		int ir = FullUpdate();
		if (!ir) return -1;
		if (ir == 2) {// immediate ua  
			BF64 w; w.bf.u64 = 0;
			for (int i = 0; i < 5; i++) {
				w |= FD[i] - zh2gxn.fsolw[i];
			}
			zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
			//cout << Char2Xout(w.bf.u64) << " seen " 
			//	<<_popcnt64(w.bf.u64 )<< endl;
			return 1;
		}
		//cout << "after update" << endl;
		//ImageCandidats();
		uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
		return (int)cc;;
	}
int ZH2_5::GoZ5G3(int fl, int* gx) {
	if (_popcnt32(fl) != 5) {
		cout << "bug fl not 5 digits" << endl;
		return -1;// not valid fl
	}
	GoZ5A(fl);// start 
	// init pm using revised gangster
	for (int i = 0; i < 54; i++) {
		int xi = C_To128[i];
		uint64_t bit = (uint64_t)1 << xi;
		if (!(cells_unsolved.bf.u64 & bit))continue;
		for (int idig = 0; idig < 5; idig++) {
			uint32_t dbit = 1 << zh2gxn.maptodigit[idig];
			if (gx[C_col[i]] & dbit) FD[idig].bf.u64 |= bit;
		}
	}

	rows_unsolved.bf.u64 = 0777777;//3*6 bits
	//ImageCandidats();

	int ir = FullUpdate();
	if (!ir) return -1;
	if (ir == 2) {// immediate gua
		BF64 w; w.bf.u64 = 0;
		for (int i = 0; i < 5; i++) {
			w |= FD[i] - zh2gxn.fsolw[i];
		}
		register uint64_t uan = w.bf.u64;
		if (zh2gxn.nkguas) {// likely redundant
			for (uint32_t i = 0; i < zh2gxn.nkguas; i++)
				if (!(~zh2gxn.knownuas[i] & (uan))) return 1;
		}
		zh2gxn.tua[zh2gxn.nua++] = uan;
		return 1;
	}
	//cout << "after update" << endl;
	////ImageCandidats();
	uint64_t cc = _popcnt64(cells_unsolved.bf.u64);
	return (int)cc;
}

int ZH2_5::DoZ5Go(int debug) { 
	//diag = debug;
	zh2gxn.nua = 0;
	//if (debug) 		ImageCandidats();
	
	/*
		if (*zh2gxn.nknownuas) {
			int ir = GetNetUaCell();
			if (ir < 0)return -1; // dead branch
			zh2b_g.guess_xcell = ir; // dead branch
		}
		*/
	ComputeNext();
	return 0;
}
inline void ZH2_5::Assign(int rdigit, int cell, int xcell) {
		FD[rdigit] &= AssignMask_Digit[cell].u64[0];
		cells_unsolved.Clear(xcell);
		int ddig = 6 * rdigit;
		rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	}
int ZH2_5::Seta(int rdigit, int xcell) { // single in cell
		int cell = From_128_To_81[xcell];
		if (FD[rdigit].Off(xcell)) return 1; // not valid
		Assign(rdigit, cell, xcell);
		BF64 *Fd = &FD[rdigit];
		BF64 * RF = &FD[4];
		for (; RF >= FD; RF--)RF->Clear(xcell);
		Fd->Set(xcell); // restore bit for digit assigned
		return 0;
	}
int ZH2_5::FullUpdate() {
		while (1) {
			if (!Update()) return 0; // game locked in update
			if (!Unsolved_Count()) return 2;
			if (ApplySingleOrEmptyCells())	return 0; // locked 
			if (zh2b_g.single_applied) 			continue;
			break;
		}
		return 1;
	}
int ZH2_5::ApplySingleOrEmptyCells() {
		zh2b_g.single_applied = 0;
		uint64_t * map = &FD[0].bf.u64;
		uint64_t unsolved = cells_unsolved.bf.u64;
		register uint64_t R2 = map[0] & map[1],	R1 = (map[0] | map[1]), 
			Map = map[2],
			R3 = R2 & Map, R4;// digits 12
		R2 |= R1 & Map; R1 |= Map;

		Map = map[3]; R4 = R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
		Map = map[4];  R4 |= R3 & Map;	R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;

		if (unsolved & (~R1)) return 1; // locked
		R1 &= ~R2;
		R1 &= unsolved; // these are new singles	
		if (R1) {
			zh2b_g.single_applied = 1;
			while (R1) {// usually a very small number of cells to assign
				uint32_t res;
				if (!bitscanforward64(res, R1)) break;
				uint64_t bit = (uint64_t)1 << res; // switch to the bit value
				R1 &= ~bit;  // clear the bit
				for (int idig = 0; idig < 5; idig++) {
					if (map[idig] & bit) {// this is the digit
						Seta(idig, res);
						goto nextr1;// finished for that cell
					}
				}
				return 1; //conflict with a previous cell assugn
			nextr1: {}
			}
			return 0;
		}
		else {
			R2 &= ~R3;
			R3 &= ~R4;
			if (!R2) R2 = R3;
			if (!R2) R2 = R4;
			// setup ok cells
			uint64_t sok = 0; // zh2gxn.unsolved_field;
			for (int i = 0; i < 5; i++) {
				sok |= (FD[i] & zh2gxn.fsolw[i]).bf.u64;
			}
			uint64_t sok2 = sok & ~cells_unsolved.bf.u64;
			sok &= cells_unsolved.bf.u64;
			int nu = *zh2gxn.nknownuas;
			//if(diag)cout << Char2Xout(sok) << " sok next cell" << endl;
			//if(diag)cout << Char2Xout(sok2) << " sok2 next cell" << endl;
				// first non hit ua to solve
			for (int iu = 0; iu < nu; iu++) {
				register uint64_t U = zh2gxn.knownuas[iu];
				if (U & sok2) continue;
				//if(diag)cout << Char2Xout(U) << "U to fill" << endl;
				U &= sok;
				if (!U)return -1; // dead branch
				if (U & R2) U &= R2;
				else if (U & R3) U &= R3;
				else U &= R4;
				bitscanforward64(zh2b_g.guess_xcell, U);
				return 0;
			}
			// all uas solved, first cell in rx
			bitscanforward64(zh2b_g.guess_xcell, R2);
			return 0;
		}
	}

void ZH2_5::Guess() {// 
		uint32_t xcell = zh2b_g.guess_xcell, cell;
		cell = From_128_To_81[xcell];
		//if (diag)			cout << "next xcell/cell="
		//	<< xcell <<" "<<cell << endl;
		uint64_t bit = (uint64_t)1 << xcell;
		int digit = zh2gxn.digit_map[zh2gxn.g0[cell]];
		if (FD[digit].bf.u64 & bit) {
			ZH2_5 * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->Seta(digit, xcell);
			mynext->ComputeNext();
		}
		for (int idig = 0; idig < 5; idig++) {
			if (FD[idig].bf.u64 & bit) {
				if (idig == digit) continue;
				ZH2_5 * mynext = this + 1; // start next guess
				*mynext = *this;
				mynext->Seta(idig, xcell);
				mynext->ComputeNext();
			}
		}
	}
void ZH2_5::ComputeNext() {
		//ImageCandidats();
		int ir = FullUpdate();
		//if (diag) {
			//cout << "compnext ir=" << ir << endl;
			//ImageCandidats();
		//}
		if (ir == 1)Guess();
		else if (ir == 2) {// solved 
			BF64 w; w.bf.u64 = 0;
			for (int i = 0; i < 5; i++) {
				w |= FD[i] - zh2gxn.fsolw[i];
			}
			uint64_t cc = _popcnt64(w.bf.u64);
			if (w.bf.u64) {
				// check no ua false
				int nu = *zh2gxn.nknownuas;
				if (nu) {
					uint64_t sok = 0; // zh2gxn.unsolved_field;
					for (int i = 0; i < 5; i++) {
						sok |= FD[i].bf.u64 & zh2gxn.fsolw[i].bf.u64;
					}
					for (int iu = 0; iu < nu; iu++) {
						register uint64_t U = zh2gxn.knownuas[iu];
						if (!(U &sok)) return;;
					}
				}
				zh2gxn.tua[zh2gxn.nua++] = w.bf.u64;
				zh2gxn.knownuas[(*zh2gxn.nknownuas)++]
					= w.bf.u64;
			}
		}
	}
int ZH2_5::Update() {
		int Shrink = 1;
		register int S, A;
		register unsigned int cl, *wcl = FD[0].bf.u32;
		while (Shrink) {
			Shrink = 0;
			if (!rows_unsolved.bf.u32[0])break;

			{register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
			if (!(AR & 077))goto digit1;

			//=digit 0
			if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
			UPDN(0, 0)	if ((AR & 7) != S) {
				AR &= 07777777770 | S;	UPWCL5(0, 2, 4, 6, 8)
			}

		digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
			UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
				AR &= 07777777707 | (S << 3);	UPWCL5(1, 3, 5, 7, 9)
			}

		digit1:	if (!(AR & 07700))goto digit2;

			if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
			UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
				AR &= 07777777077 | (S << 6); UPWCL5(0, 0, 4, 6, 8)
			}

		digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
			UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
				AR &= 07777770777 | (S << 9); UPWCL5(1, 1, 5, 7, 9)
			}

		digit2:	if (!(AR & 0770000))goto digit3;

			if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
			UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
				AR &= 07777707777 | (S << 12);	UPWCL5(0, 0, 2, 6, 8)
			}

		digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
			UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
				AR &= 07777077777 | (S << 15);	UPWCL5(1, 1, 3, 7, 9)
			}

		digit3: if (!(AR & 077000000))goto digit4;

			if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
			UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
				AR &= 07770777777 | (S << 18); UPWCL5(0, 0, 2, 4, 8)
			}

		digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
			UPDN(3, 1)if (((AR >> 21) & 7) != S) {
				AR &= 07707777777 | (S << 21); UPWCL5(1, 1, 3, 5, 9)
			}

		digit4:if (!(AR & 07700000000))goto end01234;

			if (FD[4].bf.u32[0] == CompFD[4].bf.u32[0])goto digit4b;
			UPDN(4, 0)if (((AR >> 24) & 7) != S) {
				AR &= 07077777777 | (S << 24);  UPWCL5(0, 0, 2, 4, 6)
			}

		digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
			UPDN(4, 1)if (((AR >> 27) & 7) != S) {
				AR &= 0777777777 | (S << 27);  UPWCL5(1, 1, 3, 5, 7)
			}

		end01234: rows_unsolved.bf.u32[0] = AR;
			}// end of validity for AR 01234

		}// end while

		return 1;
	}

void ZH2_5::ImageCandidats() {
		BF64  R2 = FD[0] & FD[1], R1 = FD[0] | FD[1];
		BF64 R3 = R2 & FD[2];
		R2 |= R1 & FD[2];	R1 |= FD[2];
		BF64 R4 = R3 & FD[3];
		R3 |= R2 & FD[3]; R2 |= R1 & FD[3];	R1 |= FD[3];
		BF64 R5 = R4 & FD[4];
		R4 |= R3 & FD[4]; R3 |= R2 & FD[4];
		R2 |= R1 & FD[4];	R1 |= FD[4];

		for (int i = 0; i < 6; i++) { // rows
			if ((i == 3)) {
				for (int ix = 0; ix < 54; ix++)       cout << (char)'-';
				cout << endl;
			}
			for (int j = 0; j < 9; j++) {
				if ((j == 3) || (j == 6))cout << "|";
				int cell = 9 * i + j, xcell = C_To128[cell];
				uint64_t	bit = (uint64_t)1 << xcell;
				if (!(R1.bf.u64&bit)) cout <<   "-     ";
				else if (R5.bf.u64&bit) cout << "12345 ";
				else {
					for (int i = 0; i < 5; i++)
						if (FD[i].bf.u64&bit)cout << i + 1;
					if (R4.bf.u64&bit)cout << "  ";
					else if (R3.bf.u64&bit)cout << "   ";
					else if (R2.bf.u64&bit)cout << "    ";
					else                   cout << "     ";
				}

			} // end for j
			cout << endl;
		}
		cout << endl;

	}


