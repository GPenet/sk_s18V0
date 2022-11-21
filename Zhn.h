
/*
ZhouSolver.h
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.

this is a remorphing of the code to work in recursive mode and in a 128 bit field
funtions have been added to work in generation mode
*/

/* the BIT_SET_30 pattern is defined in a band for a digit
   it is a map of candidates plus 3 bits for unsolved rows
   bits 0-8 9-17 18-26 for the rows 27-29 for unsolved rows
*/
//#include "t_128GP.h"
// tables specific to the brute force zh_tables
const extern int TblRowUnsolved[8];
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
const extern T128 zhoustart[19];
//const extern T128 AssignMask_OtherDigits[81];
 /*ZH_1D class to solve one digit 3 bands
all valid solutions  are stored in table
the table is supplied by the caller
this is somehow a private class for the brute force
and this is part of the critical code in the brute force
except for easiest puzzles immediatly solved
*/
struct ZH_1D_GLOBAL {
	BF128 *tsolw, t3; // belong to the caller
	int nsolw,lim;
	ZH_1D_GLOBAL() { t3.SetAll_1(); t3.bf.u32[3] = 0; }
	inline void Add(BF128 & s) {
		*tsolw++ = s & t3; // clear last 32 bits
		nsolw++;
	}
	int Go(BF128 & fde, BF128 *tsol,int limit=5000);
};

struct ZHOU;
/* ZH_GLOBAL2  and ZH_GLOBAL are both "static variables for the brute force
to have a better cache effect, ZH_GLOBAL is strictly limited to what is required 
for the basic brute force; 
data used in other brute force uses are in ZH_GLOBAL2
*/
struct ZH_GLOBAL2 {
	BF128 Digit_cell_Assigned[9];// init sequence
	BF128 digit_sol[9]; // final solution per digit original sort sequence
	BF128  cells_assigned;// to avoid redundancy in new assignments 
	uint64_t  npuz;
	int	xcell_to_guess;// , isfalse_on;
	GINT16 tgiven[81];
	int ngiven, digitsbf;// digitsbf to check minimum 8 digits
	int s17_b3_mini;// 17 search mode, 1 if minirows b3 not tested
	int * grid0; // using a process with known solution grid
	char * zsol,
		stdfirstsol[82],
		zerobased_sol[81];
	char  *puzfinal, *pat;
	char puz[82]; // the solved puzzle (after morph)
	void Debug();


};
struct ZH_GLOBAL { // global variables for the core brute force

	int nsol, lim, modevalid,
		icount, ntsol, single_applied,// new_single_in_Update,
		go_back,
		rdigit, loop, diag, modeguess , maxindex;
	BF128 init_3x, init_digit, pairs,triplets, pairs_naked;


	ZH_GLOBAL();
	inline void Init(int maxsols,int mvalid=0){
		nsol = go_back=0;
		lim = maxsols;
		modevalid = mvalid;
	}
	int Go_InitSudoku(char * ze);
	void ValidPuzzle(ZHOU * z);

};
/* 2 BF 128 per digit
	equivalent to F and Fcomp in Zhou
	Last 32 bits in FD[0] is a  bits field for unknown rows
	Last 32 bits in FD[1]] contain  digit mapped
*/
// class encapsulating the brute force 

#define ISFALSEON misc.bf.u32[0]
struct ZHOU{// size 32 bytes 
	BF128 FD[9][2];
	BF128 cells_unsolved,misc;// misc used to check false
//________________________________________
	int CheckValidityQuick(char *puzzle);
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	void GuessGo(int dig, BF128& s);
	void GuessInCell();
	void GuessFullDigit();
	int GuessHiddenBivalue();
	//inline void Copy(ZHOU & o);
	void Assign(int digit, int cell, int xcell);
	int Update();
	int InitSudoku(GINT16 * t, int n);
	int InitSudoku(char * zpuz);
	char * SetKnown(char * zs);
	void SetaCom(int digit, int cell, int xcell);
	inline void SetFloor(int cell){ FD[0][0] &= AssignMask_Digit[cell]; }
	inline void Seta_c(int digit, int cell){ SetaCom(digit, cell, C_To128[cell]); }
	inline void SetG(GINT16 x){ Seta_c(x.u8[1], x.u8[0]); }
	void Setcell(int cell);

	inline int Unsolved_Count(){ return cells_unsolved.Count(); }
	// standard calls
	inline void ComputeNext(){ 	if (FullUpdate())Guess();}


	// other calls and functions in Zh4_calls_variants
	inline int StatusValid(GINT16 * t, int n){
		InitSudoku(t, n); return Isvalid();
	}


	int PartialInitSudoku(GINT16 * t, int n);
	int EndInitSudoku( GINT16 * t, int n);
	//int EndInitNextUa(ZHOU & o, int bf);// 17 search check know small uas in bloc

	int IsMinimale(GINT16 * t, int n);
	//void PatFinal();
	//int GetSolvedDigitForCell(int cell);

	// inline small functions
	inline int IsOffCandidate_c(int dig, int cell){return FD[dig][0].Off_c(cell); }
	inline int IsOnCandidate_c(int dig, int cell){ return FD[dig][0].On_c(cell); }
	inline void ClearCandidate_x(int dig, int xcell){ FD[dig][0].clearBit(xcell); }
	inline void ClearCandidate_c(int dig, int cell){ FD[dig][0].clearBit(C_To128[cell]); }
	int Isvalid();
	inline int CountDigs(){
		int n =0;
		for (int i = 0; i < 9; i++)
			n += FD[i][0].Count96();
		return n;
	}

	// debugging code or print code
	void Debug(int all = 0);
	void DebugDigit(int digit);
	int GetAllDigits(int cell);
	void ImageCandidats();
	void ImageCandidats_b3();

	//==== special final check 7 search

	int PartialInitSearch17(uint32_t * t, int n);// 17 search mode
	int CallCheckB3( uint32_t bf,int nogo=0);// 17 search mode
	int CallMultipleB3x(ZHOU & o, uint32_t bf, int diag = 0);// 17 search mode
	int Apply17SingleOrEmptyCellsB3();
	int Apply17SingleOrEmptyCellsB12();
	int Full17Update();
	void Guess17(int index);
	void Compute17Next(int index) ;


 };
 
struct ZHGXN {
	BF128 fsol[9], fsolw[5];
	BF128 tua[1000];
	uint32_t floors, nua, cell_to_guess, digit_map[9];
	int * g0;

	void SetupFsol(int * grid0);
};
// class encapsulating the brute force for ua generation 2 digits
struct ZHOU2 {
	BF128 FD[2][2];
	BF128 cells_unsolved;

	int GoZ2(int fl);
	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	void Assign(int digit, int cell, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};
// class encapsulating the brute force for ua generation 2 digits
struct ZHOU3 {
	BF128 FD[3][2];
	BF128 cells_unsolved;

	int GoZ3(int fl);
	int DoZ3(int * t, int nt);
	//________________________________________
	int FullUpdate();
	int ApplySingleOrEmptyCells();
	void Guess();
	void Assign(int digit, int cell, int xcell);
	int Update();
	inline int Unsolved_Count() { return cells_unsolved.Count(); }
	void ComputeNext();
	void ImageCandidats();

};
