
#define _CRT_SECURE_NO_DEPRECATE
#define SEARCH17SOL

//#define HAVEKNOWN 
//#define DEBUGKNOWN 0  forget use sgo.vx[1]



/* program organisation
	main is the standard frame including the basic brute force 
*/

#define GTEST17_ON 1
#define UALIMSIZE 20
#define GUALIMSIZE 18
#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
//============================================== 

#include <sys/timeb.h>
#include "main.h"  // main and main tables and basic brute force
#include "go_17sol_tables.h"     
#include "Zh1b2b.h"  // brute force 2 bands  
extern SGO sgo;
//_________________ brute force handling 1 2 3 bands 
extern ZHOU    zhou[50],zhou_i;// , zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern ZH_GLOBAL2 zh_g2;
extern ZHGXN zhgxn;
extern ZHOU2 zhou2[5];
extern ZHOU3 zhou3[10];
extern ZH2B_GLOBAL   zh2b_g;
//extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands
extern ZH2GXN zh2gxn;
extern ZH2_2  zh2_2[5];
extern ZH2_3  zh2_3[10];
extern ZH2_4  zh2_4[20];
extern ZH2_5  zh2_5[20];
extern ZH2B zh2b[40], zh2b_i, zh2b_i1;
//extern ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row

FINPUT finput;
ofstream  fout1, fout2;

#include "sk_s18h.h"   //main classes of the project
#include "go_17sol_tables.h"
G17B g17b;
GEN_BANDES_12 genb12;
STD_B416 myband1, myband2;


//tables of potential bands 1+2
BF128 tbelow7[10], tbelow8[100], tbelow9[500], tbelow10[2000], tbelow11[10000];
BF128 t_7_9[2000];
uint32_t ntbelow[6], nt_7_9;//7,8,9,10,11,full clues
uint64_t tandbelow[6],tand_7_9; // same as count
uint64_t tfull[100000];

uint64_t p_cptg[40], p_cpt1g[20], p_cpt2g[100];
uint64_t p_cpt[40], p_cpt1[20];



#include "go_17_bands_cpp.h"  

#include "go_17_genb12_cpp.h"     
#include "go_17sol_bs_cpp.h"    
#include "go_17sol_commands_cpp.h"



void Go_0() {
	// open  outputs files 1.txt
	if (sgo.foutput_name) {
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_file1.txt");
		fout1.open(zn);
	}
	if (sgo.command >= 10
		&& sgo.command <20 ) {// input file expected
		if (!sgo.finput_name) {
			cerr << "missing input file name" << sgo.finput_name << endl; return;
		}
		finput.open(sgo.finput_name);
		if (!finput.is_open()) {
			cerr << "error open file " << sgo.finput_name << endl;
			return;
		}
	}
	cerr << "running command " << sgo.command << endl;
	switch (sgo.command) {
	case 0: Go_c17_00(); break; // search one band1
	case 10: Go_c17_10(); break; // run known
	case 12: Go_c17_12(); break; // band analysis and more 
	}
	cerr << "go_0 return" << endl;
}

