
/*
ZhouSolver.h
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/

/* the BIT_SET_30 pattern is defined in a band for a digit
   it is a map of candidates plus 3 bits for unsolved rows
   bits 0-8 9-17 18-26 for the rows 27-29 for unsolved rows
*/
//#include "t_128GP.h"
// tables specific to the brute force located in zh4_tables
const extern int TblRowMask[8];// rows where single  found  000 to 111
const extern int  Tblstartblock[27]; // (i%3) * 27 in zhou brute force
const extern int TblShrinkMask[512];// existing minirows 000 to 111
const extern int TblComplexMask[512]; // keep mini rows still valid optimised process
const extern int TblMaskSingle[512]; // kill in other blocks locked column /box
const extern int TblMaskDouble[512];// kill for locked in box / column
const extern int TblColumnSingle[512]; // single in column applied to shrinked bloc
const extern int TblShrinkSingle[512]; // keep only rows with single
const extern int TblRowUniq[512]; // 1 is row not defined in block  mode  to 111
const extern T128 AssignMask_Digit[81];
extern uint64_t zh2b_start[20];


/* class encapsulating the brute force with one known band
remaining clues are given in a 0-53 "cells" space
the 2 bands are located in a 64 bits field
unknown rows per digit are using 2x32 bits
first 5x6 bits for digits 0-4
second 4x6bits for digits 5-8

*/
struct ZH2B_GLOBAL { // global variables for the game table
	uint64_t * digsols; // pointer to solution grid per digit
	uint64_t ua_ret, tb1245[100], ntb1245;
	BF64 val_init1_81;// , pairs, triplets;
	BF64 Digit_cell_Assigned_init[9];
	BF64 Digit_cell_Assigned_step[9];
	BF64 fd_sols[2][9];//start puzzle/ solution
	BF64 fd_revised[9];// gangster revision of the solution
	BF64 fdsw[3][9];//morphed digits puzzle/ solution rev
	// ==============bands sols handling

	uint32_t  nua, ndigits,guess_xcell,test4b,ntest4b;
	int nsol, lim, limadd,icount, ntsol, single_applied, new_single_in_Update,
		rdigit, nctlg, go_back,diag;
	// band UA collection active band pointers and UA table to build
	int modeguess;
	int  puz0[54], gangster[9];
	int nsolw;
	char * zsol, out54[55], puzc[55];// bands 12 in output mode
	char zdebug[82];

	ZH2B_GLOBAL();
	
	void Init_g0(int* g0);
	void InitGangster(int * g0, int * g1);//common ot both GUA2s GUA3s
	uint64_t BuildUaret(BF64 * wsol);


};
/* 2 BF 128 per digit
	equivalent to F and Fcomp in Zhou
	Last 32 bits in FD[0] is a  bits field for unknown rows
	Last 32 bits in FD[1]] contain  digit mapped
*/
// class encapsulating the brute force 
struct ZH2B {// size 32 bytes 
	BF64 FD[9], CompFD[9], cells_unsolved, rows_unsolved;
	//uint64_t	cells_solved_false, active_cells,ua;
	void InitBands12(int * g0);

	//_________ brute force with uas
	void InitB1245();
	void InitB1346();
	void InitB2356();
	void InitBf(uint64_t bf);
	void InitBfG2(uint64_t bf, int c1, int d1, int c2, int d2);
	int Do4bGo();
	int Dob(int lim = 20);
	int FullUpdateb();
	int ApplySingleOrEmptyCellsb();
	inline void ComputeNextb();
	void Guessb();

	//___uas builder process


	void Init_std_bands(); // after getbands in zh2b_g
	void Init_gang(); // after getbands in zh2b_g
	inline void Assign(int digit, int cell, int xcell);
	inline int Unsolved_Count() { return rows_unsolved.Count(); }
	
	void InitTclues(uint32_t * tclues, int n);
	//void Init_2digits_banda(BF64  cellsbf);
	uint64_t IsValid(uint32_t * tclues, int n,int onlyone=0); 
	uint64_t IsValid(uint64_t bf54, int onlyone = 0);

	int Update();
	int FullUpdate();

	inline void ComputeNext();
	void GuessValidB12();
	//void GuessGo(int dig, BF64 & wsol);

	int ApplySingleOrEmptyCells();
	char * SetKnown(char * zs);
	int Seta(int digit, int xcell);
	inline int SetaC(int digit, int cell) { return Seta(digit, C_To128[cell]); }
	inline void ClearCandidate(int dig, int xcell) { FD[dig].Clear(xcell); }
	int GetAllDigits(int cell);
	void Debug(int all=0);
	void ImageCandidats();
};


struct ZH2GXN {
	uint64_t fsol[9];
	BF64 fsolw[5];
	uint64_t tua[200], unsolved_field,	uamin,uaw;

	// expand using known uas
	uint64_t  rx,r2,r3,r4,r5;// from apply single or empty cells
	uint64_t * knownuas;
	uint32_t * nknownuas,nkguas;

	uint32_t floors, nua, cell_to_guess,onlyone,
		digit_map[9], maptodigit[5], gangsters[9];
	int * g0;
	void SetupFsol(int * grid0);
	inline void InitKnown(uint64_t * kn , uint32_t *nkn) {
		knownuas = kn; nknownuas = nkn;
	}
	inline uint64_t Getsol(uint32_t bfd) {
		register uint64_t B = 0;
		for (int i = 0; i < 9; i++)if (bfd & (1 << i))
			B |= fsol[i];
		return B;
	}
};


struct ZH2_2 {
	BF64 FD[2], CompFD[2], cells_unsolved,
		rows_unsolved;
	uint64_t diag;
	void GoZ2A(int fl);
	int GoZ2G2(int fl, int c1, int d1, int c2, int d2);
	void DoZ2Go();
	int GoZ2D(int  fl);

	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	inline void Assign(int digit, int cell, int xcell);
	int Seta(int digit, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};

struct ZH2_3 {
	BF64 FD[3], CompFD[3], cells_unsolved,
		rows_unsolved;
	uint64_t diag;
	void GoZ3A(int fl);
	int GoZ3(int fl);
	int GoZ3G2(int fl,int c1,int d1,int c2,int d2);
	int DoZ3Go(int debug=0);
	int GoZ3G3(int fl, int * gx, int debug = 0);// see main process

	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	inline void Assign(int digit, int cell, int xcell);
	int Seta(int digit, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};
struct ZH2_4 {
	BF64 FD[4], CompFD[4], cells_unsolved, rows_unsolved;
	void GoZ4A(int fl);
	int GoZ4(int fl);
	int GoZ4G2(int fl, int c1, int d1, int c2, int d2);
	int GoZ4G3(int fl, int* gx);// see main process
	int DoZ4Go();
	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	inline void Assign(int digit, int cell, int xcell);
	int Seta(int digit, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};
struct ZH2_5 {
	BF64 FD[5], CompFD[5], cells_unsolved, rows_unsolved;
	//uint64_t diag;

	void GoZ5A(int fl);
	int GoZ5(int fl);
	int GoZ5G2(int fl, int c1, int d1, int c2, int d2);
	int GoZ5G3(int fl, int* gx);// see main process
	int DoZ5Go(int debug=0);
	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	int GetNetUaCell();
	void Guess();
	inline void Assign(int digit, int cell, int xcell);
	int Seta(int digit, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};





