// GPengine.h: interface for the CGPengine class.
//
//////////////////////////////////////////////////////////////////////
//W.B.Langdon @ cs.ucl.ac.uk 23 August 2000 Elvis Hand-Eye cordination experiment
//Changes
//WBL 11 Jul 2025 Add entropy
//WBL 25 Jul 2025 Reduce MaxInstr to 1e5
//WBL  4 Jul 2025 Add Simplify InstrLen2 Instr2 Changed
//WBL 29 Jun 2025 for stats add output
//WBL 28 Jun 2025 For long runs fitness as int, not array remove FitnessCase->Output(), testval
//WBL 19 Jun 2025 For 2004 Memorial Uni Mackey-Glass experiments
//WBL 19 Jun 2025 For Linux g++ 11.5.0
//WBL 23 Aug 00 Add input file specification class, change FitnessCase
//              More registers
//              Change DefaultPopSize, MaxInstr, IntRangeXxx

#if !defined(AFX_GPENGINE_H__2330E6E0_6CEE_11D3_B1A6_0060086C3065__INCLUDED_)
#define AFX_GPENGINE_H__2330E6E0_6CEE_11D3_B1A6_0060086C3065__INCLUDED_

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <fstream>
#include <stdio.h>
//#include <conio.h>
#include <assert.h>

using namespace std;

#define BOOL int

extern int nthreads;

#ifdef elvis
const int GenerateLimit = 125000; //1000000;
const int DefaultPopSize = 5000;
const int MaxInstr  = 50; //25;
const int PCrossover = 90;
const int PMutation = 80;

typedef double retval;
typedef double OP;
#else
/*Mackey Glass chaotic time series geccolb.tex r1.20 gp.details*/
//const int GenerateLimit = 500*(100000-1);
const int GenerateLimit = 500*(500-1);
const int DefaultPopSize = 500;
const int MaxInstr  = 4000000;
const int MaxInitInstr = 14;
const int PCrossover = 90;
const int PMutation = 40;

typedef unsigned char testval;
typedef unsigned char retval;
typedef int OP;
#endif

typedef OP instr[4];
typedef char var[2];

#define validvar(x)    (x>= 0 && x < NumVar)
#define validreg(x)    (x > IntRangeEnd && x <= (IntRangeEnd + NumVar))
#define validopcode(x) (x>= 0 && x < NumOp)
#ifdef elvis
#define validvalue(x)  1
#else
#define validvalue(x)  (x >= IntRangeStart && x <= IntRangeEnd)
#endif

#define validInstr(x) (\
  validvar(   x[0]) && \
  validreg(   x[1]) && \
  validopcode(x[2]) && \
  (validreg(x[3]) || validvalue(x[3])))

struct Individual
{
	int InstrLen;
	int InstrLen2;
	int fitness;
	instr *Instr;
	instr *Instr2;
	retval *output; //for convergence stats
//	double* fitness;
};


class __FitnessCase
{
	testval* data;
	int records;
	int size;
	__FitnessCase(const int r, const int n): records(r), size(n) {
		data = new testval[r*size];}
	~__FitnessCase() { delete[] data;}
	int Num() const  { return records; }
	int Size() const { return size; } 
	void Data(const int r, const int x, const testval d) {
		assert(0<=r && r<records && 0<=x && x<size);
		data[r*size+x] = d;
	}
	testval Data(const int r, const int x) const {
		assert(0<=r && r<records && 0<=x && x<size);
		return data[r*size+x];
	}
	friend class _FitnessCase;
};
class _FitnessCase
{
	__FitnessCase* input;
	__FitnessCase* output;
public:
	_FitnessCase(const int r, const int in, const int out) {
		input  = new __FitnessCase(r,in);
		output = new __FitnessCase(r,out);
	}
	~_FitnessCase() { delete output; delete input; }
	int Num() const { return input->Num(); }
	void    Input(const int r, const int x, const testval d) { input->Data(r,x,d); }
	testval  Input(const int r, const int x) const    { return input->Data(r,x);   }
	void   Output(const int r, const int x, const testval d){ output->Data(r,x,d); }
	testval Output(const int r) const                { return output->Data(r,0);   }
	int     Input() const { return input->Size(); }
	void write(const char* filename) const;
};

// Constants available for a individual
#ifdef elvis
const int IntRangeStart = -100;
const int IntRangeEnd = 100;
#else
const int IntRangeStart = 0;
const int IntRangeEnd = 127;
#endif

const int NumVar = 8;

/*
const int a = 0+IntRangeEnd+1;
const int b = 1+IntRangeEnd+1;
const int c = 2+IntRangeEnd+1;
const int d = 3+IntRangeEnd+1;


const int NumOp = 4;

#define plus 0 
#define minus 1
#define mul 2
#define div 3
*/

class dataspec {
	int size;
	int max;
	struct chain {chain* next; int pos; };
	chain* first;
	//chain* last = null;
public:
	dataspec(): size(0), max(-1), first(NULL) {;};
	BOOL add(const int pos);
	~dataspec();
	int Pos(const int i) const;
	int Size() const;
	int Max() const;
};

class CGPengine  
{
#ifdef elvis
	enum { plus, minus, mul, div, _acos, _asin, _atan2, _cos, _fabs, _hypot, _sign, _sin, _sqrt, _tan, NumOp };
#else
	enum { plus, minus, mul, div, NumOp };
#endif

public:

	CGPengine();
	CGPengine(int PopSize);
	virtual ~CGPengine();
	
	int Init();
   void LoadTrainingData(const char *filename, const dataspec& input, const dataspec& output);
	int GeneratePop();
    int displayDiffEval(ostream& out, const retval* evals1, const int instrlen1, const int instr1, const retval* evals2, const int instrlen2, const int instr2, const int display) const;
   void displayEval(ostream& out,
		  /*const*/ Individual &parent1, const int start1, const int len1,
		  /*const*/ Individual &parent2, const int start2, const int len2,
		      const int mutations[4],
		  /*const*/ Individual &child) const;
   void Evolve(ostream&);
   void GenerateCode(ostream&) const;
   void GenerateBestCode(const int/*not in use,ostream&*/) const;
   void GenerateBestCode(ostream& out) const;
   void write(const char* filename) const {FitnessCase->write(filename);};
   void RunFitness(ofstream& out) const;
   void LoadRun();
  //private:
	
	// tournament selection
       void Tournament(); //const int objective was used by elvis
	
	int CrossLenOk(const int InstrLen1, const int xolen1, const int xolen2) const;
	int ChooseXO(const int InstrLen, int& start) const;
       void ChooseXO(const int InstrLen1, int& start1, int& len1,
		     const int InstrLen2, int& start2, int& len2) const;
       void CrossPart(const int Winner, const int start1, const int len1, const int start2, const int Looser);
       void Cross1(const int Winner1, const int start1, const int len1, const int Winner2, const int start2, const int len2, const int Looser2);
       void Crossover();
       void Mutation();
	int MutateInstr(int Instr,Individual &I);
       void Reproduction();

   void CalcFitness(Individual &I);
   void display(const int test, const retval reg[NumVar]) const;
   void display(const int test, const retval reg[NumVar], const int line, const instr &ii) const;
    int Changed(const int looser, const int winner);
    int change2(const int looser, const int winner) const; //Instr2    differ
    int changeg(const int looser, const int winner) const; //genotypes differ
   void Simplify(Individual &I, int* Needed = NULL) /*const if Needed given*/;
   void Interpret(const Individual &I, retval* output2) const;
   void Interpret(const int InstrLen, const instr *Instr, const int fcase, retval &output, retval* output2 = NULL) const;

   void GenerateIndividual(Individual &I);
    int InstrArg(const OP code) const;
 retval InstrArg(const OP code, const retval reg[]) const;
    int InstrReg(const OP code) const;
 retval InstrReg(const OP code, const retval reg[]) const;
     OP GenerateReg() const;
     OP GenerateArg() const;
   void GenerateInstr(instr I);
		
	void GenerateCode(const instr &i, ostream&) const;
	void GenerateCode(const Individual &I, const int* needed /*=NULL*/, ostream&) const;
	void GenerateHead(const int, ostream&) const;
	void GenerateTail(ostream&) const;
	char *Var(int op) const;
	void Var(const OP i, ostream& out) const;
	void VarOrVal(const OP i, ostream& out) const;
	//char *Int(int i);
	char *Op(const int i) const;

	Individual *Pop;
	int PopSize;
	var *VarList;
	
	_FitnessCase *FitnessCase;
	//int NumFitnessCases;
	
	// tournament index of individs
	int Winner1,Winner2,Looser1,Looser2;

	int BestFitness;
//	int	   BestIndividual;
	
};

int  SizeTrainingData(const char *filename,const int columns);
BOOL TrainingData(FILE* infile, const int columns, double row[]);

void   clear();
void   increment(const uint64_t key);
double entropy();
  

#endif // !defined(AFX_GPENGINE_H__2330E6E0_6CEE_11D3_B1A6_0060086C3065__INCLUDED_)
