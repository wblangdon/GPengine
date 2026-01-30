#if defined __skylake_avx512__ 
extern const char* rev = "$Revision: 1.129 $ AVX512";
#include <immintrin.h>
#else
extern const char* rev = "$Revision: 1.129 $";
#endif
#if defined __skylake__
error not written yet
#endif
// GPengine.cpp: implementation of the CGPengine class.
//
//////////////////////////////////////////////////////////////////////
//W.B.Langdon @ cs.ucl.ac.uk 23 August 2000 Elvis Hand-Eye cordination experiment
//Changes
//WBL 14 Dec 2025 apply simplified avx_20251118_1763483201.dif
//WBL 11 Dec 2025 Add __skylake_avx512__ Interpret64 TODO add Magpie patch
//WBL 23 Sep 2025 Add rev
//WBL 13 Aug 2025 Ensure get full entropy trace (stats.awk limits file size)
//WBL 12 Aug 2025 Rerun mg.code0 again but with entropy and only 1st example
//    unify displayEval's out with cout and filter to reduce volume of output
//WBL  8 Aug 2025 Expand trace to include bad children in last generation
//WBL  7 Aug 2025 trace disruption in last generation
//WBL  1 Aug 2025 Want to trace FDP, start with "stats" to record genetic and fitness changes
//WBL 25 Jul 2025 Try using perf for elapsed time
//WBL 24 Jul 2025 instrument Simplify for gpfunc
//WBL 22 Jul 2025 Add pthreads. Overhead 1st 300 gens 0.37microseconds per Interpret call
//WBL 21 Jul 2025 Revert to r1.41, ie put AVX to one side (r1.65)
//WBL  4 Jul 2025 Add Simplify InstrLen2 Instr2 Changed evals
//WBL 28 Jun 2025 for stats add xolong, output
//WBL 28 Jun 2025 For long runs remove displayfitness(out.. Fitness as int, testval
//WBL 25 Jun 2025 add gpfunc conditional code
//WBL 24 Jun 2025 bugfix Best fitness calculation
//WBL 21 Jun 2025 bugfix Reproduction
//WBL 19 Jun 2025 For 2004 Memorial Uni Mackey-Glass experiments
//WBL 19 Jun 2025 For Linux g++ 11.5.0
//WBL 23 Aug 00 TODO change popsize, 
//                   
//              read from file, Add more variables, Add more variables, 
//              non-linear opcodes (eg logics, trig), 


/* compile:
  g++ -fpermissive -fmax-errors=2 -c GPengine.cpp 
 */

#include "GPengine.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include <climits>
#include <sys/timeb.h>
#include <pthread.h>
#ifdef gpfunc
//Based on https://github.com/wblangdon/linux_perf_api/blob/main/demo_perf.c
//and test_prog.c r1.4
#include <unistd.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
static long
perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
		int cpu, int group_fd, unsigned long flags){
  int ret;
  ret = syscall(__NR_perf_event_open, hw_event, pid, cpu,group_fd, flags);
  return ret;
}
int perf_start(){
  struct perf_event_attr pe;
  int fd;

  memset(&pe, 0, sizeof(pe));
  //set up perf's event data structure.
  //There are many things perf can report, here ask for instruction count
  //Note your CPU may not support all perf's options
  //pe.type = PERF_TYPE_HARDWARE;	//PERF_TYPE_HW_CACHE;
  pe.type = PERF_TYPE_SOFTWARE;
  pe.size = sizeof(pe);
  //PERF_COUNT_HW_CACHE_MAX
  //PERF_COUNT_HW_CACHE_OP_MAX
  //PERF_COUNT_HW_CACHE_RESULT_MAX
  //? not yet ? PERF_COUNT_SW_MAX
  //pe.config = PERF_COUNT_HW_INSTRUCTIONS; //PERF_COUNT_HW_CACHE_MISSES;
  pe.config = PERF_COUNT_SW_CPU_CLOCK;
  pe.disabled = 1;
  pe.exclude_kernel = 1;
  pe.exclude_hv = 1;

  fd = perf_event_open(&pe, 0, -1, -1, 0);
  if (fd == -1) {
    fprintf(stderr, "Error on perf_event_open %llx\n", pe.config);
    return EXIT_FAILURE;
  }

  ioctl(fd, PERF_EVENT_IOC_RESET, 0);
  ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
  return fd;
}
void perf_end(const int fd, const char* text) {
  long long count = -1;
  ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
  if(errno != 0) {//check all is well 
    //Might get for example errno == EBADF 9 /* Bad file number */ bad fd ?
    fprintf(stderr,"%s Opps errno is %d\n",text,errno);
    exit(errno);
  }
  read(fd, &count, sizeof(count)); //instructions is first 8 bytes read from fd
  close(fd);

  printf("%s took %lld ticks\n",text,count);
}

#endif /*gpfunc*/

unsigned int rrand(int max) 
{
	return (rand()%max);	
}

int randtest(double percent) 
{
	return (unsigned) int(10.0*percent)>rrand(1000u);
}

int xolong = 0; //just for reporting stats
int evals  = 0; //just for reporting stats

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CGPengine::CGPengine()
{	
	PopSize = DefaultPopSize;
	FitnessCase = NULL;
	Pop = NULL;
}

CGPengine::CGPengine(int NewSize)
{
	PopSize = NewSize;		
	FitnessCase = NULL;
	Pop = NULL;
}

//from magpie/avx512 eval_perf.cpp r1.28
uint32_t div32[256*256];
void init_div(){
  for(int op1=0;op1<256;op1++){
  for(int op2=0;op2<256;op2++){
    if (op2 == 0)
      div32[op1*256+op2] = 0;
    else
      div32[op1*256+op2] = (retval)(op1/op2);
  }}
}
//end  eval_perf.cpp

int CGPengine::Init()
{
	int Error=0;

	// create variable list
	VarList = new var[NumVar];
	
	for (int i=0;i<NumVar;i++)
	{
		VarList[i][0] = char(int('a')+i);
		VarList[i][1] = '\0';
	}
#ifdef stats
	BestFitness = INT_MAX;
#endif //*stats*/
//	BestIndividual = -1;

	init_div();

	return Error;
}


//WBL re-write to check for formating errors
int  SizeTrainingData(const char *filename,const int columns)
{
	FILE* infile = NULL;
	infile = fopen(filename,"r");
	if(infile != NULL) 	{printf("reading from %s ", filename); fflush(stderr); }
	else {
		fprintf(stderr, "fopen r %s failed\n", filename);
		exit(1);
	}
	BOOL eof=false;
	int line = 0;
	int count = 0;
	while(eof==false) {
		line++;
		char buf[513]; 
		int c;
		int i;for(i=0;(i<sizeof(buf)-1)&&((c=getc(infile))!=EOF)&&(c!=10);i++) buf[i]=char(c);
		buf[i]=0;
		eof=(c==EOF);
		//printf("%s %d\n",buf,eof);
		assert(columns<=16);
		float data[16];
		const int n = sscanf(buf,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
			&data[0],&data[1],&data[2], &data[3], &data[4], &data[5], &data[6], &data[7],
			&data[8],&data[9],&data[10],&data[11],&data[12],&data[13],&data[14],&data[15]);
		if(n==0) {}//assume comment ignore
		else if(n==EOF && eof) {} //end loop on end of file
        else if(n<columns || n>16) {
			fprintf(stderr,"Bad data line %d file %s %d '%s' \n",line,filename,n,buf);
		} 
		else count++;
	}//endloop
    fclose(infile);
	printf(" %d records \n", count);
	return count;
}
//------------------------------------------------------------------------------------
BOOL TrainingData(FILE* infile, const int columns, double row[])
{
	BOOL eof=false;
	int line = 0;
	while(eof==false) {
		line++;
		char buf[513]; 
		int c;
		int i;for(i=0;(i<sizeof(buf)-1)&&((c=getc(infile))!=EOF)&&(c!=10);i++) buf[i]=char(c);
		buf[i]=0;
		eof=(c==EOF);
		//printf("%s %d\n",buf,eof);
		assert(columns<=16);
		float data[16];
		const int n = sscanf(buf,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
			&data[0],&data[1],&data[2], &data[3], &data[4], &data[5], &data[6], &data[7],
			&data[8],&data[9],&data[10],&data[11],&data[12],&data[13],&data[14],&data[15]);
		if(n==0) {}//assume comment ignore
        else if(n<columns || n>16) {
			fprintf(stderr,"Bad data line %d input file '%s' \n",line,buf);
		} 
		else {
			for(i=0;i<columns;i++) row[i] = double(data[i]);
			return true;
		}
	}//endloop
	return false;
}
//------------------------------------------------------------------------------------
void CGPengine::LoadTrainingData(const char *filename, const dataspec& input, const dataspec& output)
{
	if(input.Size()>NumVar||output.Size()>NumVar)
		fprintf(stderr,"Warning: #inputs %d #outputs %d exceeds #registers %d \n",input.Size(),output.Size(),NumVar);
	if(output.Size() != 1) {
		fprintf(stderr,"Error only one output now supported %d\n",output.Size());
		exit(1);
	}

	const int max = (input.Max() > output.Max())? input.Max() : output.Max(); 
	const int num = SizeTrainingData(filename,max);
	FitnessCase = new _FitnessCase(num,input.Size(),output.Size());

	FILE* infile = fopen(filename,"r");
	double* row = new double[max];
	for (int i=0;i<FitnessCase->Num();i++)
	{
		if(!TrainingData(infile,max,row)) exit(1);
		for(int j=0; j<input.Size(); j++) {
		  const int     p = input.Pos(j)-1;
		  const testval d = row[p];
		  if(d != row[p]) {
		    cerr<<"LoadTrainingData("<<filename<<",,) i="<<i<<" j="<<j<<" "
			<<"precision error on input row["<<p<<"] "<<row[p]<<endl;
		    exit(1);
		  }
		  FitnessCase-> Input(i,j,d);
		}
		for(int j=0; j<output.Size();j++) {
		  const int     p = output.Pos(j)-1;
		  const testval d = row[p];
		  if(d != row[p]) {
		    cerr<<"LoadTrainingData("<<filename<<",,) i="<<i<<" j="<<j<<" "
			<<"precision error on output row["<<p<<"] "<<row[p]<<endl;
		    exit(1);
		  }
		  FitnessCase->Output(i,j,d);
		}
	}
	delete[] row;
	fclose(infile);
	
/*	ifstream file(filename);
	
	float inp,outp;

	file >> NumFitnessCases;

	FitnessCase = new _FitnessCase[NumFitnessCases];

	for (int i=0;i<NumFitnessCases;i++)
	{
			file >> inp;
			file >> outp;
			FitnessCase[i].input = double(inp);
			FitnessCase[i].output = double(outp);
	}

	file.close();

	return Error;
*/
}

CGPengine::~CGPengine()
{

	if (FitnessCase) 
		delete FitnessCase;


	if (Pop)
	{
		for (int i=0;i<PopSize;i++)
			delete [] Pop[i].Instr;
	
		delete Pop;
	}
}


int CGPengine::GeneratePop()
{
	int Error=0;

	
	(time(NULL)%65535);
	
	// create population
	Pop = new Individual[PopSize];

	for (int i=0;i<PopSize;i++)
		GenerateIndividual(Pop[i]);
	
	return Error;
}

/*
/////////////////////////////////////////////////////////////////////////////////////////////////////
// heap.cc   test heapsort 275--284 Foundations of Computer Science, C Edition, Aho and Ullman, 1995
// W.B.Langdon cwi.nl 26 Sep 1999
//#define main_version "Revision: 1.00 "

//Modifications (in reverse order)

//WBL 26 Sep 99  

//#include "pch.h"


void swap(int A[], const int i, const int j) {
  const int temp = A[i];
  A[i] = A[j];
  A[j] = temp;
}

void bubbleDown(int A[], const int i, const int n) {
  int child = 2*i;
  if(child <  n && A[child+1] > A[child]) child++;
  if(child <= n && A[i]       < A[child]) {
    swap(A,i,child);
    bubbleDown(A,child,n);
  }
}

void heapify(int A[], const int n) {
  for(int i = n/2; i>=1; i--) bubbleDown(A,i,n);
}

void deletemax(int A[], int *pn) {
  swap(A, 1,*pn);
  --(*pn);
  bubbleDown(A,1,*pn);
}

void heapsort(int A[], const int n) {
  heapify(A,n);
  int i = n;
  while(i>1) deletemax(A,&i);
}

#define MAX 100
int A[MAX+1];
int main(int argc, char * argv[])
{
  int n,x;
  n=0;
  while(n<MAX && scanf("%d",&x)!=EOF) A[++n]=x;

  heapsort(A,n);

  for(int i=1;i<=n;i++) cout<<A[i]<<endl;

  return 1;
}
*/
//end heap.cc /////////////////////////////////////////////////////////////////////////////////////////
#ifndef elvis
/*
//from sym.cc
// W. B. Langdon cs.bham.ac.uk Revision: 1.25 
// This code is released for non-commercial use only
// W.Langdon cs.bham.ac.uk 

//Modifications (reverse order):
// WBL 25 Apr 2000  For ROVL
#endif
class sort {
private:
  retval data; int y;
public:
  inline int Line() const {return y;};
  inline void set(const retval a, const int b)
    { data=a; y=b; }
  inline BOOL operator<(const sort& b) const 
    { return data < b.data; };
  inline BOOL operator>(const sort& b) const 
    { return data > b.data; };
  friend ostream& operator << (ostream& os, const sort& t);
};
ostream& operator << (ostream& os, const sort& t) {
  os<<t.data<<" "<<t.y; 
  return os;
}


inline void swap(sort A[], const int i, const int j) {
  const sort temp = A[i];
  A[i] = A[j];
  A[j] = temp;
}

inline void bubbleDown(sort A[], const int i, const int n) {
  int child = 2*i;
  if(child <  n && A[child+1] > A[child]) child++;
  if(child <= n && A[i]       < A[child]) {
    swap(A,i,child);
    bubbleDown(A,child,n);
  }
}

inline void heapify(sort A[], const int n) {
  for(int i = n/2; i>=1; i--) bubbleDown(A,i,n);
}

inline void deletemax(sort A[], int *pn) {
  swap(A, 1,*pn);
  --(*pn);
  bubbleDown(A,1,*pn);
}

inline void heapsort(sort A[], const int n, const int top) {
  heapify(A,n);
  int i = n;
  while(/* i> 1***  i>top) deletemax(A,&i);
}
#ifndef elvis
cf r1.29 */
#endif /*elvis*/
//end sym.cc /////////////////////////////////////////////////////
int FitnessCaseNum = 0;
#ifdef stats
//last crossovers
int start1 = -1;
int start2 = -1;
int len1   = -1;
int len2   = -1;
int mutations[8]; //last mutations
int fdp_example = 0;

size_t fcn_nv = 0;
inline retval eval(const size_t t, const size_t r, const retval* evals, const size_t instrlen, const size_t instr) {
  assert(t>=0 && t<FitnessCaseNum);
  assert(r>=0 && r<NumVar);
  const size_t i = t*instrlen*NumVar+instr*NumVar+r;
  assert(i>=0 && i < instrlen*fcn_nv);
  return evals[i];
}
  
#define Evals1(t,r) eval(t,r,evals1,instrlen1,instr1)
#define Evals2(t,r) eval(t,r,evals2,instrlen2,instr2)
inline int nodiff_eval(const retval* evals1, const int instrlen1, const int instr1, const retval* evals2, const int instrlen2, const int instr2){//sanity check only
  for(int i=0;i<FitnessCaseNum;i++){
    for(int r=0;r<NumVar;r++){
      if(Evals1(i,r) != Evals2(i,r)) return 0;
    }
  }
  return 1;
}
double eval_entropy(const retval* evals, const size_t instrlen, const size_t instr) {
  clear();
  for(int i=0;i<FitnessCaseNum;i++){
    uint64_t reg = 0;
    for(int r=0;r<NumVar;r++){
      const uint64_t e  = eval(i,r,evals,instrlen,instr);
      const uint64_t ee = e << r*8;
      reg |= ee;
    }
    increment(reg);
  }
  return entropy();
}
int CGPengine::displayDiffEval(ostream& out, const retval* evals1, const int instrlen1, const int instr1, const retval* evals2, const int instrlen2, const int instr2, const int display) const{
  int diff0=0; //Only R0
  int diffs=0;
  int diffa=0; //count all eval differences
  for(int i=0;i<FitnessCaseNum;i++){
    if(Evals1(i,0) != Evals2(i,0)) diff0++;
    for(int r=0;r<NumVar;r++){
      if(Evals1(i,r) != Evals2(i,r)) {diffs++; break;}
    }
    for(int r=0;r<NumVar;r++){
      if(Evals1(i,r) != Evals2(i,r)) diffa++;
    }
  }
  char buff[80];
  sprintf(buff,"%4d ",diff0);out<<buff;
  sprintf(buff,"%4d ",diffs);out<<buff;
  sprintf(buff,"%4d ",diffa);out<<buff;
  sprintf(buff,"%10.7f ",eval_entropy(evals1,instrlen1,instr1));out<<buff;
  sprintf(buff,"%10.7f ",eval_entropy(evals2,instrlen2,instr2));out<<buff;
  if(fdp_example > 1) return diffa;//long form only once
  if(diffs==0 && display==0) return diffa;
  out<<" ";
  int d=0;
  for(int i=0; d<2 && i<FitnessCaseNum;i++){
    int dd = (diffs==0); //display forces output or look for differences
    for(int r=0; dd==0 && r<NumVar;r++){
      if(Evals1(i,r) != Evals2(i,r)) dd=1;
    }
    if(dd){
      d++;
      char buff[7];sprintf(buff,"T%-4d",i);out<<buff;
      for(int r=0;r<NumVar;r++){
	const int k = i*NumVar+r;
	if(Evals1(i,r) != Evals2(i,r)) {
	  char buff[10];sprintf(buff,"%3d!=%-3d ",Evals1(i,r),Evals2(i,r));out<<buff;
	}
	else{
	  char buff[10];sprintf(buff,"%3d      ", Evals1(i,r));out<<buff;
	}
      }
    }
  }
  out<<flush;
  return diffa;
#undef Evals2
#undef Evals1
}
size_t output2_size = 0; //only for debug
//follow argument order in cross1 Looser2
void CGPengine::displayEval(ostream& out,
		  /*const*/ Individual &parent1, const int start1, const int len1,
		  /*const*/ Individual &parent2, const int start2, const int len2,
		      const int mutations[4],
	          /*const*/ Individual &child) const {
  /*
  cout<<"displayEval(out,"
      <<"parent1("<<parent1.InstrLen<<" "<<parent1.InstrLen2<<", "<<parent1.fitness<<"),"
      <<start1<<","<<len1<<" "
      <<"parent2("<<parent2.InstrLen<<" "<<parent2.InstrLen2<<", "<<parent2.fitness<<"),"
      <<start2<<","<<len2<<" "
      <<"mutations[4],"
      <<"child("  <<child.InstrLen<<" "  <<child.InstrLen2  <<", "<<child.fitness  <<")"
      <<endl;
  */
  out <<"displayEval(out,"
      <<"parent1("<<parent1.InstrLen<<" "<<parent1.InstrLen2<<", "<<parent1.fitness<<"),"
      <<start1<<","<<len1<<" "
      <<"parent2("<<parent2.InstrLen<<" "<<parent2.InstrLen2<<", "<<parent2.fitness<<"),"
      <<start2<<","<<len2<<" "
      <<"mutations[4],"
      <<"child("  <<child.InstrLen<<" "  <<child.InstrLen2  <<", "<<child.fitness  <<")"
      <<endl;

  fcn_nv = FitnessCaseNum*NumVar;

  const size_t evals1_size = ((size_t)parent1.InstrLen)*fcn_nv;
#ifdef NDEBUG
  retval* evals1 = (retval*)malloc(evals1_size*sizeof(retval));
#else
  retval* evals1 = (retval*)calloc(evals1_size,sizeof(retval));
#endif /*NDEBUG*/
  if(evals1) ;//out<<"ok"<<endl;
  else {out<<"evals1_size "<<evals1_size<<" failed"<<endl; return;}
  output2_size = evals1_size;
  Interpret(parent1, evals1);
  output2_size = 0;
  int* needed1 = (int*)calloc(parent1.InstrLen,sizeof(int));
  Simplify(parent1,needed1);
  int cneeded1=0;for(int i=0;i<parent1.InstrLen;i++){if(needed1[i])cneeded1++;}
  if(cneeded1 != parent1.InstrLen2){
    out<<"ERROR cneeded1 "<<cneeded1<<" != parent1.InstrLen2 "<<parent1.InstrLen2
       <<endl;}
  
  const size_t evals2_size = ((size_t)parent2.InstrLen)*fcn_nv;
  //char buf[80];
  //sprintf(buf,"evals2_size %lld %zx parent2.InstrLen %d",evals2_size,evals2_size,parent2.InstrLen);
  //out<<buf<<endl;
  retval* evals2 = (retval*)calloc(evals2_size,sizeof(retval));
  if(evals2); //out<<"ok"<<endl;
  else {out<<"evals2_size "<<evals2_size<<" failed"<<endl;
    free(needed1);
    free(evals1);
    return;}
  output2_size = evals2_size;
  Interpret(parent2, evals2);
  output2_size = 0;

  const size_t evals_size = ((size_t)child.InstrLen)*fcn_nv;
  retval* evals = (retval*)calloc(evals_size,sizeof(retval));
  if(evals) ;//out<<"ok"<<endl;
  else {out<<"evals_size "<<evals_size<<" failed"<<endl;
    free(evals2);
    free(needed1);
    free(evals1);
    return;}

  output2_size = evals_size;
  Interpret(child,   evals);
  output2_size = 0;
  int* needed  = (int*)calloc(child.InstrLen, sizeof(int));
  Simplify(child,needed);
  int cneeded=0;for(int i=0;i<child.InstrLen;i++){if(needed[i])cneeded++;}
  if(cneeded != child.InstrLen2){
    out<<"ERROR cneeded "<<cneeded<<" != child.InstrLen2 "<<child.InstrLen2
       <<endl;}
  /*if(fdp_example == 1){//long form only once
  out<<"\nWinner2 "<<endl;
  GenerateCode(parent2,NULL,out);
  out<<"\nWinner1 "<<endl;
  GenerateCode(parent1,needed1,out);
  out<<endl;
  }*/
  int* mut = (int*)calloc(child.InstrLen,sizeof(int));
  for(int i=0;i<4;i++){
    const int k = mutations[i];
    if(k >= 0){
      if(!(k>=0 && k<child.InstrLen)){
	cerr<<"k "<<k<<" out of bounds "<<child.InstrLen<<" child.InstrLen\n";
	exit(99);
      }
      mut[k] = 1;
    }
  }
  int write = 0;
  for(int i=0;i<child.InstrLen;i++){
    if(fdp_example >= 1) write = 1; //to get entropy always use long form
    if(mut[i])           write = 1;
    if(i==start2)        write = 1;
    if(i==start2+len1)   write = 1;
    int p = i;
    int q = -1;
    if(i<start2) {}
    else if(i<start2+len1) { q = i-start2+start1; assert(q>=0 && q < parent1.InstrLen);} //P2
    else { p = i-len1+len2; assert(p>=0 && p < parent2.InstrLen);} //P1
    if(write == 0) {
      if(q == -1){
	assert(nodiff_eval(evals2,parent2.InstrLen,p,evals,child.InstrLen,i));}//sanity check only
      else {
	assert(nodiff_eval(evals1,parent1.InstrLen,q,evals,child.InstrLen,i));}//sanity check only
    } //end if check no output is needed
    if(write == 0) continue;
    char buff[11];sprintf(buff,"%7d ",i);out<<buff;
    GenerateCode(child.Instr[i],out); out<<" ";
    if(i<start2) {
      out<<"P1 "; GenerateCode(parent2.Instr[i],out);}
    else if(i<start2+len1){
      out<<"P2 "; GenerateCode(parent1.Instr[q],out);}
    else { p = i-len1+len2; assert(p>=0 && p < parent2.InstrLen);
      out<<"P1 "; GenerateCode(parent2.Instr[p],out);}

    out<<" ";
    out<<((mut[i])?     "M" : "_"); const int display = (i+1<child.InstrLen)? mut[i+1] : 0;
    out<<((needed1[i])? "N" : "_");
    out<<((needed[i])?  "n" : "_");
    if(q == -1){
      out<<" ";if(!displayDiffEval(out,evals2,parent2.InstrLen,p,evals,child.InstrLen,i,display)) write = 0;
    } else {
      out<<" ";if(!displayDiffEval(out,evals1,parent1.InstrLen,q,evals,child.InstrLen,i,display)) write = 0;
    }
    out<<endl;
    if(fdp_example < 1 && write == 0) out<<endl;
  }
  out<<endl;

  free(mut);
  free(needed);
  free(evals);
  free(evals2);
  free(needed1);
  free(evals1);
}
#endif /*stats*/

void displayfitness(ostream & out,
					const int index, const Individual& I, 
					const double mean, const int max	) {
	out<<index<<" "<<" "<<I.fitness;
	out<<" #Instr "<<I.InstrLen;
	out<<" "/*<<min[j]<<" "*/<<mean<<" "<<max;
}

/* too much even for NDEBUG
void sanity(const Individual &I, const _FitnessCase* FitnessCase){//Test if fitness and output match
  int fitness = 0;
  for (int i=0;i<FitnessCase->Num();i++){
    const int error = FitnessCase->Output(i) - I.output[i];
    fitness += error * error; //should not overflow 32bit int
  }
  if(I.fitness != fitness) {
    fprintf(stderr,"sanity output %p I.fitness %d != %d\n",I.output,I.fitness,fitness);
    exit(99);
  } else {
    fprintf(stderr,"sanity output %p I.fitness %d ok\n",I.output,I.fitness);
  }
  fflush(NULL);
  assert(I.fitness == fitness);
  }
*/
void CGPengine::GenerateBestCode(const int cnt/*not in use, ostream& out*/) const
{
/* too much even for NDEBUG
	{for (int i=0;i<PopSize;i++) {
	    printf("Pop[%3d] %6d %7d\n: ",i,Pop[i].InstrLen,Pop[i].fitness);
	  for(int j=0;j<FitnessCase->Num();j++) printf("%3d ",Pop[i].output[j]);
	  printf("\n");
	  sanity(Pop[i],FitnessCase);
	  if(i>0) assert(Pop[i].output != Pop[0].output); //something broken somewhere
	}}
	cout <<flush;
*/
//WBL 29 Aug 2000 rewrite to sort population

	cout << cnt << " Fitness: "<<flush;
/*make rendering easier for emacs will need fixup if define elvis
#ifdef elvis
	June 2025 following appears problemantic when FitnessCase->Output() == 1

	sort** scores = new sort*[FitnessCase->Output()];
	{for(int j=0;j<FitnessCase->Output();j++) scores[j] = new sort[PopSize];}
	int**  ranks  = new int*[FitnessCase->Output()];
	{for(int j=0;j<FitnessCase->Output();j++){
	   ranks[j]  = new int[PopSize];
#ifndef NDEBUG
	   memset(ranks[j],0xff,sizeof(int)*PopSize);
#endif
	}}
	
	double* sum  = new double[FitnessCase->Output()]; {for(int j=0;j<FitnessCase->Output();j++) sum[j] = 0;}
	double* min  = new double[FitnessCase->Output()]; {for(int j=0;j<FitnessCase->Output();j++) min[j] = FLT_MAX;}
	double* max  = new double[FitnessCase->Output()]; {for(int j=0;j<FitnessCase->Output();j++) max[j] = -FLT_MAX;}
	double  ssum = 0; 
	double  smin = FLT_MAX;
	double  smax = -FLT_MAX;

	for (int i=0;i<PopSize;i++) {
		double sf = 0;
		for(int j=0;j<FitnessCase->Output();j++) {
			const double f = Pop[i].fitness[j];
			scores[j][i].set(f,i);
			sum[j] += f;
			if(f<min[j]) min[j] = f;
			if(f>max[j]) max[j] = f;
			sf     += f;
		}
		sf /= FitnessCase->Output();
		ssum += sf;
		if(sf<smin) smin = sf;
		if(sf>smax) smax = sf;
	}//end scan whole population
	//const int top = PopSize; //for timebeing
	{for(int j=0;j<FitnessCase->Output();j++) {
		heapsort(&scores[j][-1], PopSize, 1); 
	}}
	{for(int j=0;j<FitnessCase->Output();j++) {
		for (int i=0;i<PopSize;i++) {
#ifndef NDEBUG
			assert(ranks[j][scores[j][i].Line()] == -1);
#endif
			ranks[j][scores[j][i].Line()] = i;
		}//end scan whole population
	}}

	int BestRank = INT_MAX;
	int Best = 0;
	{for (int i=0;i<PopSize;i++) {
		int rank = 0;
		for(int j=0;j<FitnessCase->Output();j++) {
#ifndef NDEBUG
			assert(ranks[j][scores[j][i].Line()] >= 0);
#endif
			rank += ranks[j][i];
		}
		if(rank<BestRank) {
			BestRank = rank;
			Best = i;
		}
	}}//Overal best is one with lowest sum of ranks;

	if(FitnessCase->Output()>1) {
		assert(Best >= 0 && Best < PopSize);
		const double srank = double(BestRank)/FitnessCase->Output();
		displayfitness(cout,FitnessCase->Output(),Best,Pop[Best],ssum/PopSize,smax); cout<<" "<<srank<<endl;
		out<<"//";
		displayfitness( out,FitnessCase->Output(),Best,Pop[Best],ssum/PopSize,smax);  out<<" "<<srank<<endl;
		GenerateHead(cnt,out);
		GenerateCode(Pop[Best],out);
		GenerateTail(out);
		out<<endl;
	}

	{for(int j=0;j<FitnessCase->Output();j++) {
		const int first = scores[j][0].Line();
		const Individual Best = Pop[first];
		const double mean = sum[j]/PopSize;
		cout<< "                ";
		displayfitness(cout,FitnessCase->Output(),first,Best,mean,max[j]); cout<<endl;

		out<<"//"<<j<<" ";
		displayfitness(out, FitnessCase->Output(),first,Best,mean,max[j]); out<<endl;
		GenerateHead(cnt,out);
		GenerateCode(Best,out);
		GenerateTail(out);
	}}
#endif /*elvis*/

//	if(FitnessCase->Output()!=1){
//	  cerr<<"FitnessCase->Output()="<<FitnessCase->Output()<<" not 1\n";
//	  exit(99);
//	}

	int   min  =  INT_MAX;
	int   max  =  INT_MIN;
	int long long sum = 0;
	int   lmin =  INT_MAX;
	int   lmax =  INT_MIN;
	long long int lsum = 0;
	int  lmin2 =  INT_MAX;
	int  lmax2 =  INT_MIN;
	long long int lsum2 = 0;
	int Best   = -1;
	int numBest = 0;
	{for (int i=0;i<PopSize;i++) {
	    const int f = Pop[i].fitness;
	    sum += f;
	    if(f < min){ min = f; Best = i; numBest = 0;}
	    if(f== min)  numBest++;
	    if(f > max)  max = f;
	    const int l = Pop[i].InstrLen;
	    lsum += l;
	    if(l < lmin) lmin = l;
	    if(l > lmax) lmax = l;
	    const int l2 = Pop[i].InstrLen2;
	    lsum2 += l2;
	    if(l2 < lmin2) lmin2 = l2;
	    if(l2 > lmax2) lmax2 = l2;
	}}
	assert(Best >= 0);
	displayfitness(cout,Best,Pop[Best],sum/PopSize,max);
	cout<<" "<<    lmin<<" "<<((double)lsum) /PopSize<<" "<<lmax<<flush;
	cout<<"\t"<<   lmin2<<" "<<((double)lsum2)/PopSize<<" "<<lmax2<<flush;
	int equiv = 0;
	{for (int i=0;i<PopSize;i++) {
	 for (int j=0;j<i;j++) {
	   const int e = (memcmp(Pop[i].output,Pop[j].output,FitnessCase->Num()*sizeof(retval)) == 0);
	   if(e) equiv++;
	}}}
	int eBest = 0;
	{for (int j=0;j<PopSize;j++) {
	   const int e = (memcmp(Pop[Best].output,Pop[j].output,FitnessCase->Num()*sizeof(retval)) == 0);
	   if(e) eBest++;
	}}
	cout<<",\t"<<numBest<<" "<<eBest<<" "<<equiv<<" "<<xolong;
	cout<<" "<<evals;
	cout<<endl;
	assert(numBest > 0);
	assert(  eBest > 0);
	assert(  eBest <= numBest);
	xolong = 0;
	evals  = 0;
/*
	out<<"//"<<" ";
	displayfitness( out,1,Best,Pop[Best],sum/PopSize,max);
	out <<" "<<    lmin<<" "<<((double)lsum)/PopSize<<" "<<lmax;
	out <<endl;
	GenerateHead(cnt,out);
	GenerateCode(Pop[Best],out);
	GenerateTail(out);

/*make rendering easier for emacs will need fixup if define elvis
#ifdef elvis
	delete[] max;
	delete[] min;
	delete[] sum;
	{for(int j=0;j<FitnessCase->Output();j++) delete[] ranks[j];}
	{for(int j=0;j<FitnessCase->Output();j++) delete[] scores[j];}
#endif /*elvis**
	cout<<flush;
	out<<endl<<endl;
*/
}
void CGPengine::GenerateBestCode(ostream& out) const
{
  int   min  =  INT_MAX;
  int   max  =  INT_MIN;
  int long long sum = 0;
  int Best   = -1; //Resolve ties conviently
  for (int i=0;i<PopSize;i++) {
    const int f = Pop[i].fitness;
    sum += f;
    if(f < min){ min = f; Best = i; }
    if(f > max)  max = f;
  }
  assert(Best >= 0);
  displayfitness(out,Best,Pop[Best],sum/PopSize,max); out<<endl;
  out<<endl;
  GenerateCode(Pop[Best],NULL,out);
}

void CGPengine::GenerateCode(ostream& out) const
{
	for (int i=0;i<PopSize;i++)
	{
		/*
		// Header
		cout << "double Individual(char v[])" << endl << "{" << endl;
		cout << "double a=b=c=d=0;" << endl << "a=v[0];"  << endl;
		*/
		GenerateCode(Pop[i],NULL,out);
		// footer
		//cout << "return a;" << endl;
	}
}



#ifdef stats
int CGPengine::changeg(const int looser, const int winner) const { //genotype change for stats only
  if(Pop[looser].InstrLen == Pop[winner].InstrLen &&
     memcmp(Pop[looser].Instr, Pop[winner].Instr, sizeof(instr)*Pop[winner].InstrLen) == 0) {
    return 0;
  }
  return 1;
}
int CGPengine::change2(const int looser, const int winner) const {//as Changed but no side effects
  if(Pop[looser].InstrLen2 == Pop[winner].InstrLen2 &&
     memcmp(Pop[looser].Instr2, Pop[winner].Instr2, sizeof(instr)*Pop[winner].InstrLen2) == 0) {
    return 0;
  }
  return 1;
}
#else /*stats*/
int CGPengine::change2(const int looser, const int winner) const {//dummy to keep linker happy
  return 0;
}
#endif /*stats*/
int CGPengine::Changed(const int looser, const int winner){
  if(Pop[looser].InstrLen2 == Pop[winner].InstrLen2 &&
     memcmp(Pop[looser].Instr2, Pop[winner].Instr2, sizeof(instr)*Pop[winner].InstrLen2) == 0) {
    Pop[looser].fitness = Pop[winner].fitness;
    memcpy(Pop[looser].output,Pop[winner].output,sizeof(retval)*FitnessCase->Num());
    return 0;
  }
  return 1;
}

void CGPengine::Evolve(ostream& out)
{
	int cnt=0;
	
	while (/*!kbhit() && */ (cnt<=GenerateLimit || GenerateLimit<=0) /*&& BestFitness-Pop[BestIndividual].InstrLen>0*/)
	{
		Tournament();//cnt%FitnessCase->Output());

		// crossover??
		if (randtest(PCrossover))
		{
#ifdef stats
			cout << "cnt " << cnt <<" "<<flush;
			memset(mutations,255,8*sizeof(int));
#endif
			Crossover();		

			// mutate loosers??
			if (randtest(PMutation)) {
#ifdef stats
				cout << "cnt " << cnt <<" "<<flush;
#endif
				Mutation();		
			}

			// calulate fittness for loosers
			Simplify(Pop[Looser1]);
			const int C11 =           Changed(Looser1,Winner1);
			const int C12 = (C11==0)? change2(Looser1,Winner2) : Changed(Looser1,Winner2);
			if(C11 && C12) CalcFitness(Pop[Looser1]);
			Simplify(Pop[Looser2]);
			const int C21 =           Changed(Looser2,Winner1);
			const int C22 = (C21==0)? change2(Looser2,Winner2) : Changed(Looser2,Winner2);
			if(C21 && C22) CalcFitness(Pop[Looser2]);
#ifdef stats
			const int c11 = (memcmp(Pop[Winner1].output,Pop[Looser1].output,FitnessCase->Num()*sizeof(retval)) != 0);
			const int c12 = (memcmp(Pop[Winner1].output,Pop[Looser2].output,FitnessCase->Num()*sizeof(retval)) != 0);
			const int c21 = (memcmp(Pop[Winner2].output,Pop[Looser1].output,FitnessCase->Num()*sizeof(retval)) != 0);
			const int c22 = (memcmp(Pop[Winner2].output,Pop[Looser2].output,FitnessCase->Num()*sizeof(retval)) != 0);
			cout << "cnt " << cnt <<flush;
			cout << " W1 "<<Winner1<<" len "<<Pop[Winner1].InstrLen<<" "<<Pop[Winner1].InstrLen2<<" fit "<< Pop[Winner1].fitness<<flush;
			cout << " W2 "<<Winner2<<" len "<<Pop[Winner2].InstrLen<<" "<<Pop[Winner2].InstrLen2<<" fit "<< Pop[Winner2].fitness<<flush;
			cout << " L1 "<<Looser1<<" len "<<Pop[Looser1].InstrLen<<" "<<Pop[Looser1].InstrLen2
			     <<" diff " << changeg(Looser1,Winner1)<<" "<<changeg(Looser1,Winner2)
			     <<" diff2 "<< C11 <<" "<< C12
			     <<" fit "<< Pop[Looser1].fitness
			     <<" diffo "<<c11<<" "<<c21<<flush;
			cout << " L2 "<<Looser2<<" len "<<Pop[Looser2].InstrLen<<" "<<Pop[Looser2].InstrLen2
			     <<" diff " << changeg(Looser2,Winner1)<<" "<<changeg(Looser2,Winner2)
			     <<" diff2 "<< C21 <<" "<< C22
			     <<" fit "<< Pop[Looser2].fitness
			     <<" diffo c12 "<<c12<<" c22 "<<c22<<flush;
			cout << endl;

			//all Looser2 crossover children in last generation
			if(cnt > GenerateLimit-PopSize /*&& fdp_example == 0*/){
		      //if(Pop[Looser2].fitness <= BestFitness){
		      //if(C22 && c22==0){
			  fdp_example++;
			  cout<<"Example FDP "<<fdp_example<<" "
			      <<Winner2<<"(len "<<Pop[Winner2].InstrLen<<" "<<Pop[Winner2].InstrLen2<<" fit "<< Pop[Winner2].fitness<<" ) "
			      <<Winner1<<"(len "<<Pop[Winner1].InstrLen<<" "<<Pop[Winner1].InstrLen2<<" fit "<< Pop[Winner1].fitness<<" ) "
			      <<" child "
			      <<Looser2<<"(len "<<Pop[Looser2].InstrLen<<" "<<Pop[Looser2].InstrLen2<<" fit "<< Pop[Looser2].fitness<<" ) "
			      <<"start1="<<start1<<" "
			      <<"start2="<<start2<<" "
			      <<"len1="<<len1<<" "
			      <<"len2="<<len2<<" "
			      <<"C22 "<<C22<<" "<<"c22 "<<c22<<" "
			      <<"mutations[ "<<flush;
			  for (int i=0;i<4;i++) cout<< mutations[i+4] << " " ;
			  cout<<"]"<<endl;
			  /*
			  out <<"Example FDP "<<fdp_example<<" "
			      <<Winner2<<"(len "<<Pop[Winner2].InstrLen<<" "<<Pop[Winner2].InstrLen2<<" fit "<< Pop[Winner2].fitness<<" ) "
			      <<Winner1<<"(len "<<Pop[Winner1].InstrLen<<" "<<Pop[Winner1].InstrLen2<<" fit "<< Pop[Winner1].fitness<<" ) "
			      <<" child "
			      <<Looser2<<"(len "<<Pop[Looser2].InstrLen<<" "<<Pop[Looser2].InstrLen2<<" fit "<< Pop[Looser2].fitness<<" ) "
			      <<"start1="<<start1<<" "
			      <<"start2="<<start2<<" "
			      <<"len1="<<len1<<" "
			      <<"len2="<<len2<<" "
			      <<"C22 "<<C22<<" "<<"c22 "<<c22<<" "
			      <<"mutations[ "<<flush;
			  for (int i=0;i<4;i++) out << mutations[i+4] << " " ;
			  out <<"]"<<endl;
			  */
			  //out<<"\nWinner2 "<<Winner2<<endl;
			  //GenerateCode(Pop[Winner2],out);
			  //out<<"\nWinner1 "<<Winner1<<endl;
			  //GenerateCode(Pop[Winner1],out);
			  //out<<endl;
			  //globals start2,len2 set by last Crossover for Looser2
			  //Pickup from last Mutation() of Looser2
			  displayEval(cout,
				      Pop[Winner1],start1,len1,Pop[Winner2],start2,len2,
				      &mutations[4],Pop[Looser2]);
			}//}}
#endif /*stats*/
		}
		// no crossover, just reproduct winners
		else
			Reproduction();
		
		cnt+=2;

		if (cnt%PopSize<2)
		{	
			if(cnt % (PopSize*100) < 2){
			  const time_t t1 = time(NULL);
			  cout << cnt <<" "<< ctime(&t1) << flush;
			}
			GenerateBestCode(cnt);
		
		}
			
	}
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//												PRIVATE														//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void CGPengine::Tournament() //const int obj)
{

	// pick four individs from pop
	int ind1,ind2,ind3,ind4;
	
	ind1 = rrand(PopSize);
	do 
	{
		ind2 = rrand(PopSize);
	} while (ind2==ind1);
				
	do
	{
		ind3 = rrand(PopSize);
	} while ((ind3==ind1)||(ind3==ind2));
				
	do 
	{
		ind4 = rrand(PopSize);
	} while ((ind4==ind1)||(ind4==ind2)||(ind4==ind3));

	// tournament
	if (Pop[ind1].fitness<Pop[ind2].fitness)
	{
		Winner1=ind1;
		Looser1=ind2;
	}
	else
	{
		Winner1=ind2;
		Looser1=ind1;
	}
	
	if (Pop[ind3].fitness<Pop[ind4].fitness)
	{
		Winner2=ind3;
		Looser2=ind4;
	}
	else
	{
		Winner2=ind4;
		Looser2=ind3;
	}
/*
	/*For simplicity assume only one output** assert(obj==0);
	if (Pop[Winner1].fitness[obj] <= BestFitness)
	{
		printf("Tournament(%d) %f BestIndividual %d updated Winner1 %d %f \n",
		       cnt,BestFitness,BestIndividual,Winner1,Pop[Winner1].fitness[obj]);
		fflush(NULL);
		BestFitness = Pop[Winner1].fitness[obj];
		BestIndividual = Winner1;
	}
	
	if (Pop[Winner2].fitness[obj] <= BestFitness)
	{
		printf("Tournament(%d) %f BestIndividual %d updated Winner2 %d %f \n",
		       cnt,BestFitness,BestIndividual,Winner2,Pop[Winner2].fitness[obj]);
		fflush(NULL);
		BestFitness = Pop[Winner2].fitness[obj];
		BestIndividual = Winner2;		
	}
	
	assert(0<=BestIndividual && BestIndividual<PopSize &&
		   BestIndividual != Looser1 && BestIndividual != Looser2);
*/
}	

void CGPengine::CrossPart(const int Winner, const int start1, const int len1, const int start2, const int Looser)
{
  assert(start1 >= 0);
  assert(start1      <= Pop[Winner].InstrLen);
  assert(len1   >= 0);
  assert(start1+len1 <= Pop[Winner].InstrLen);
  assert(start2 >= 0);
  assert(start2+len1 <= MaxInstr);

		// copy winner segment into looser
		memcpy(Pop[Looser].Instr[start2],Pop[Winner].Instr[start1],sizeof(instr)*len1);
}
void CGPengine::Cross1(const int Winner1, const int start1, const int len1, const int Winner2, const int start2, const int len2, const int Looser2)
{
  assert(start1 >= 0 && start1 < Pop[Winner1].InstrLen);
  assert(len1   >= 0 && len1  <= Pop[Winner1].InstrLen);
  assert(start2 >= 0 && start2 < Pop[Winner2].InstrLen);
  assert(len2   >= 0 && len2  <= Pop[Winner2].InstrLen);
  
		// copy winner2 to looser 2
		CrossPart(Winner2,0,start2,0,Looser2);
		
		// copy winner1 segment into looser 2		
		CrossPart(Winner1,start1,len1,start2,Looser2);

		// append rest of winner2
		const int len = Pop[Winner2].InstrLen-(start2+len2);
		CrossPart(Winner2,start2+len2,len,start2+len1,Looser2);
			
		Pop[Looser2].InstrLen = start2+len1+len;		

  assert(len >= 0);
  assert(Pop[Looser2].InstrLen > 0 && Pop[Looser2].InstrLen<=MaxInstr);
}

int CGPengine::CrossLenOk(const int InstrLen1, const int xolen1, const int xolen2) const
{
	assert(xolen1  <= InstrLen1);
	assert(xolen1  >= 0 && xolen1   <= MaxInstr);
	assert(xolen2  >= 0 && xolen2   <= MaxInstr);
	assert(InstrLen1 > 0 && InstrLen1 <= MaxInstr);

	const int newlen = InstrLen1 + xolen2 - xolen1;
	return newlen > 0 && newlen <= MaxInstr;
}
int CGPengine::ChooseXO(const int InstrLen, int& start) const
{
	const int s1 = rrand(InstrLen);
	const int s2 = rrand(InstrLen);
	if(s1<=s2) { //s1==s2 should be ok
	  start = s1;
	  return s2-s1;
	} else {
	  start = s2;
	  return s1-s2;
	}
}
void CGPengine::ChooseXO(const int InstrLen1, int& start1, int& len1,
			 const int InstrLen2, int& start2, int& len2) const
{
	len1 = ChooseXO(InstrLen1, start1);
	len2 = ChooseXO(InstrLen2, start2);
	if(CrossLenOk(InstrLen1,len1,len2) && CrossLenOk(InstrLen2,len2,len1)) return;

	//violating MaxInstr choose len1,len2 again subject to start1,start2 etc. and MaxInstr
	xolong++; //just for reporting stats
	const int max1   = InstrLen1 - start1;
	const int max2   = InstrLen2 - start2;
	const int choice = (max1 > max2);
	const int InstrLena = (choice)? InstrLen1 : InstrLen2;
	const int InstrLenb = (choice)? InstrLen2 : InstrLen1;
	const int lenb   = (choice)? len2 : len1;
	const int max    = (choice)? max1 : max2;
	const int rx     = rrand(1+max);
	int lena         = rx;
	const int childa = (InstrLena+(lenb-lena) > MaxInstr);
	if(childa) lena  = lenb + InstrLena - MaxInstr;
	const int childb = (InstrLenb+(lena-lenb) > MaxInstr);
	if(childb) lena  = lenb + MaxInstr  - InstrLenb;

	if(choice){
	  len1 = lena;
	} else {
	  len2 = lena;
	}

	/*cout<<"ChooseXO("<<InstrLen1<<","<<start1<<","<<len1<<","
                         <<InstrLen2<<","<<start2<<","<<len2<<") "
	    <<"choice="<<choice<<" "
	    <<"childb too long="<<childb<<" "
	    <<"childa too long="<<childa<<" "
	    <<"lena="<<lena<<" "
	    <<"lenb="<<lenb<<" "
	    <<"rx="<<rx<<" "
	    <<"lena="<<lena<<" "
	    <<endl;
	if(CrossLenOk(InstrLen1,len1,len2) == 0)
	  cout<<"CrossLenOk(Pop[Winner1]."<<InstrLen1<<","<<len1<<","<<len2<<") failed\n";
	if(CrossLenOk(InstrLen2,len2,len1) == 0)
	  cout<<"CrossLenOk(Pop[Winner2]."<<InstrLen2<<","<<len2<<","<<len1<<") failed\n";
	*/
	assert(len1 >= 0);
	assert(start1+len1 <= InstrLen1);
	assert(start1+len1 <= MaxInstr);
	assert(len2 >= 0);
	assert(start2+len2 <= InstrLen2);
	assert(start2+len2 <= MaxInstr);

	assert(InstrLen1+(len2-len1) <= MaxInstr);
	assert(InstrLen2+(len1-len2) <= MaxInstr);
}
void CGPengine::Crossover()
{
	// select randomplace in winners
#ifndef stats
	int start1 = -1;
	int start2 = -1;
	int len1   = -1;
	int len2   = -1;
#endif /*not stats*/
	ChooseXO(Pop[Winner1].InstrLen,start1,len1,
		 Pop[Winner2].InstrLen,start2,len2);
	assert(CrossLenOk(Pop[Winner1].InstrLen,len1,len2));
	assert(CrossLenOk(Pop[Winner2].InstrLen,len2,len1));

	// Create looser2 
	Cross1(Winner1,start1,len1,Winner2,start2,len2,Looser2);

	// Create looser 1	
	Cross1(Winner2,start2,len2,Winner1,start1,len1,Looser1);

#ifdef stats
	cout << "Crossover"<<flush;
	cout << " Winner1 "<<Winner1<<" len "<<Pop[Winner1].InstrLen<<" "<<Pop[Winner1].InstrLen2<<" fit "<< Pop[Winner1].fitness<<flush;
	cout << " Winner2 "<<Winner2<<" len "<<Pop[Winner2].InstrLen<<" "<<Pop[Winner2].InstrLen2<<" fit "<< Pop[Winner2].fitness<<flush;
	cout << " start1 "<<start1<<" len1 "<< len1<< " start2 "<< start2<<" len2 "<<len2<<flush;
	cout << " Looser1 "<<Looser1<<" len "<<Pop[Looser1].InstrLen<<flush;
	cout << " Looser2 "<<Looser2<<" len "<<Pop[Looser2].InstrLen<<flush;
	cout << endl;
#endif /*stats*/
}


void CGPengine::Mutation()
{
#ifdef stats
	cout << "Mutation "<<Looser1 << " " << Looser2 << " " <<flush;
#endif /*stats*/

	for (int i=0;i<4;i++)
	{
		// select instruction to mutate
		int Instr1 = rrand(Pop[Looser1].InstrLen);
		int Instr2 = rrand(Pop[Looser2].InstrLen);
#ifdef stats
		const int c1 = MutateInstr(Instr1,Pop[Looser1]);
		const int c2 = MutateInstr(Instr2,Pop[Looser2]);
		cout << Instr1 << " " << c1 << " ; " <<flush;
		cout << Instr2 << " " << c2 << " , " <<flush;
		assert(i+4 < 8);
		mutations[i]   = (c1)? Instr1 : -1;
		mutations[i+4] = (c2)? Instr2 : -1;
#else
		MutateInstr(Instr1,Pop[Looser1]);
		MutateInstr(Instr2,Pop[Looser2]);	
#endif /*stats*/
	}
#ifdef stats
	cout << endl;
#endif /*stats*/
}

int CGPengine::MutateInstr(int Instr,Individual &I)
{
#ifdef stats
	instr old;
	assert(sizeof(instr) == 16);
	memcpy(old,I.Instr[Instr],sizeof(instr));
#endif /*stats*/
	// Select which pos to mutate
	int pos = rrand(4);
	int val;

	switch(pos)
	{
		// change var
		case 0:
				do
				{
					val = rrand(NumVar);
				} while (val == I.Instr[Instr][0]);
				I.Instr[Instr][0] = val;
			break;
		// change var/constant
		case 1: 
				I.Instr[Instr][1] = GenerateReg();
				break;
		// change operator
		case 2:
				do
				{
					val = rrand(NumOp);
				} while (val == I.Instr[Instr][2]);
				I.Instr[Instr][2] = val;				
				break;
		// change var/constant
		case 3: 
				I.Instr[Instr][3] = GenerateArg();
				break;

		default:
			break;
	}
#ifdef stats
	return memcmp(old,I.Instr[Instr],sizeof(instr)) != 0;
#else
	return 0;
#endif /*stats*/
}

void CGPengine::Reproduction()
{
	// Looser 1
	Pop[Looser1].fitness   = Pop[Winner1].fitness;
	Pop[Looser1].InstrLen  = Pop[Winner1].InstrLen;
	Pop[Looser1].InstrLen2 = Pop[Winner1].InstrLen2;
	memcpy(Pop[Looser1].Instr, Pop[Winner1].Instr, sizeof(instr)*Pop[Winner1].InstrLen);
	memcpy(Pop[Looser1].Instr2,Pop[Winner1].Instr2,sizeof(instr)*Pop[Winner1].InstrLen2);
	memcpy(Pop[Looser1].output,Pop[Winner1].output,sizeof(retval)*FitnessCase->Num());

	// Looser 2
	Pop[Looser2].fitness   = Pop[Winner2].fitness;
	Pop[Looser2].InstrLen  = Pop[Winner2].InstrLen;
	Pop[Looser2].InstrLen2 = Pop[Winner2].InstrLen2;
	memcpy(Pop[Looser2].Instr, Pop[Winner2].Instr, sizeof(instr)*Pop[Winner2].InstrLen);
	memcpy(Pop[Looser2].Instr2,Pop[Winner2].Instr2,sizeof(instr)*Pop[Winner2].InstrLen2);
	memcpy(Pop[Looser2].output,Pop[Winner2].output,sizeof(retval)*FitnessCase->Num());
}

retval* Reg = NULL;
retval* Reg64 = NULL;
retval* Output = NULL;

#if defined __skylake__ || defined __skylake_avx512__ 
#define regtype unsigned char
#define npar 64
enum { plus_op, minus_op, mul_op, div_op, NumOp };

//from magpie/avx512 experiments eval.cpp Revision: 1.20
//WBL include file XML template for eval_perf.cpp Magpie experiments Revision: 1.20

//Modifications:
//WBL 18 Nov 2025 Remove code not used by regtype8
//WBL 13 Nov 2025 Allow uint32_t
//WBL 12 Nov 2025 Allow uint16_t
//WBL 11 Nov 2025 Use -Dregtype in manual AVX512 64-byte conversion
//WBL 13 Oct 2025 based on GPengine.cpp r1.65 Interpret16 etc

int InstrArg(const OP code) //return code not value
{
	const int x = code-IntRangeEnd-1;
	return x;
}

__m512i InstrArg16(const OP code, const uint8_t registers[], const int j)
{
	if(code>IntRangeEnd){
		const int x = InstrArg(code);
	const __m128i a    = _mm_loadu_epi8(&registers[x*npar+j]);
	const __m512i aa   = _mm512_cvtepi8_epi32(a); //sign extends
	const __m512i zero = _mm512_set1_epi32(0);
	const __mmask64 k  = 0xeeeeeeeeeeeeeeee;
	return _mm512_mask_blend_epi8(k, aa, zero);
	}
	else {
		return _mm512_set1_epi32(code);
	}
}
__m512i InstrArg32(const OP code, const uint8_t registers[], const int j)
{
	if(code>IntRangeEnd){
		const int x = InstrArg(code);
	const __m256i a    = _mm256_loadu_epi8(&registers[x*npar+j]);
	return _mm512_cvtepi8_epi16(a); //sign extends
	//const __m512i aa   = _mm512_cvtepi8_epi16(a); //sign extends
	//const __m512i zero = _mm512_set1_epi32(0);
	//const __mmask64 k  = 0;//0xaaaaaaaaaaaaaaaa;
	//return _mm512_mask_blend_epi8(k, aa, zero);
	}
	else {
		return _mm512_set1_epi16(code);
	}
}
//Sanity check: have removed code for uint16_t uint32_t
static_assert(sizeof(regtype) == sizeof(uint8_t));
__m512i InstrArg64(const OP code, const uint8_t registers[])
{
	if(code>IntRangeEnd){
		const int x = InstrArg(code);
		return _mm512_loadu_epi8(&registers[x*npar]);
	}
	else {
		return _mm512_set1_epi8(code);
	}
}
int InstrReg(const OP code) //return code not value
{
	const int x = code-IntRangeEnd-1;
	return x;
}

__m512i InstrReg16(const OP code, const uint8_t registers[], const int j)
{
	const int x = InstrReg(code);
	const __m128i a    = _mm_loadu_epi8(&registers[x*npar+j]);
	const __m512i aa   = _mm512_cvtepi8_epi32(a); //sign extends
	const __m512i zero = _mm512_set1_epi32(0);
	const __mmask64 k  = 0xeeeeeeeeeeeeeeee;
	return _mm512_mask_blend_epi8(k, aa, zero);
}
__m512i InstrReg32(const OP code, const uint8_t registers[], const int j)
{
	const int x = InstrReg(code);
	const __m256i a    = _mm256_loadu_epi8(&registers[x*npar+j]);
	return _mm512_cvtepi8_epi16(a); //sign extends
	//const __m512i aa   = _mm512_cvtepi8_epi16(a); //sign extends
	//const __m512i zero = _mm512_set1_epi32(0);
	//const __mmask64 k  = 0;//0xaaaaaaaaaaaaaaaa;
	//return _mm512_mask_blend_epi8(k, aa, zero);

}
__m512i InstrReg64(const OP code, const uint8_t registers[])
{
	const int x = InstrReg(code);
	return _mm512_loadu_epi8(&registers[x*npar]);
}

//based on CGPengine::Interpret (without stats)
void Interpret64(const int InstrLen, const instr *Instr, regtype registers[NumVar*npar], const uint32_t div32[256*256]) {
static_assert(npar==64);
// run program
for (int i=0;i<InstrLen;i++) {		
  if(Instr[i][2]>=div_op) {
  //case div:
  for(int j=0;j<npar;j+=16){
    const __m512i a          = InstrReg16(Instr[i][1],registers,j);
    const __m512i a256       = _mm512_slli_epi32(a, 8); //shift left by 1 byte
    const __m512i by         = InstrArg16(Instr[i][3],registers,j);
    const __m512i index_vec  = _mm512_add_epi32(a256,by);
    const __m512i result_vec = _mm512_i32gather_epi32(index_vec, div32, 4);
    const __m128i c          = _mm512_cvtepi32_epi8(result_vec);
    _mm_store_si128((__m128i*)&registers[Instr[i][0]*npar+j], c);
  }//endfor j=0,16,32,48
  } else {
    switch(Instr[i][2])
      {
      case plus_op:{
	const __m512i a = InstrReg64(Instr[i][1],registers);
	const __m512i b = InstrArg64(Instr[i][3],registers);
	__m512i c = _mm512_add_epi8 (a, b);
	_mm512_store_epi32((__m512i*)&registers[Instr[i][0]*npar],c);
	break;
      }
      case minus_op:{
	const __m512i a = InstrReg64(Instr[i][1],registers);
	const __m512i b = InstrArg64(Instr[i][3],registers);
	__m512i c = _mm512_sub_epi8 (a, b);
	_mm512_store_epi32((__m512i*)&registers[Instr[i][0]*npar],c);
	break;
      }
      case mul_op:{
	for(int j=0;j<npar;j+=32){
	  const __m512i a  = InstrReg32(Instr[i][1],registers,j);
	  const __m512i b  = InstrArg32(Instr[i][3],registers,j);
	  const __m512i c  = _mm512_mullo_epi16(a, b);
	  const __m256i cc = _mm512_cvtepi16_epi8(c);
	  _mm256_storeu_si256((__m256i_u*)&registers[Instr[i][0]*npar+j],cc);
	}//endfor j=0,32
	break;
      }
      }//end switch (+ - *)
  }//endif/else div
 }
}//end Interpret64
//end eval.cpp Revision: 1.20
#endif /*SSE or AVX*/

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
int slot; //protect with mutex
CGPengine*  engine; //only one instance of CGPengine in use
Individual* newguy; //only one Individual at a time
void* interpret(void* junk) {
  for (;;) {
    {const int e = pthread_mutex_lock(&mutex);   assert(e==0);}
#if defined __skylake__ || defined __skylake_avx512__ 
    const int i = slot;   slot += npar;
#else
    const int i = slot++;
#endif
    {const int e = pthread_mutex_unlock(&mutex); assert(e==0);}
    if(i >= FitnessCaseNum) break;
    const Individual* I = newguy;
#if defined __skylake__ || defined __skylake_avx512__ 
    __attribute__((aligned(sizeof(__m512i))))
    regtype registers[NumVar*npar];
    int ntests;
    {//hide I,i
	// registers
	assert(sizeof(regtype) == sizeof(retval));
	const int fcase = i;
	const int I = fcase*8;
	assert(NumVar == 8);
	assert(fcase >= 0 && fcase < 1201);
	assert(I>=0 && I+7 < 8*1201);
	assert(I>=0 && I+63 < 9728);

	memcpy(registers,&Reg64[I],NumVar*npar*sizeof(retval));
	ntests = (fcase+npar < 1201)? npar : 1201 - fcase;
#ifndef NDEBUG
	for (int fcase=i;fcase<i+npar && fcase<1201;fcase++) {
	  for (int j=0;j<NumVar;j++) {
	    const int i = fcase*NumVar + j;
	    assert(i >= 0 && i < 9728);
	    assert( ((i<128*8 && Reg[i] >=0) || Reg[i] > 27) && Reg[i] < 170);
	  }
	}//endfor fcase
#endif
    }//end scope
    // run program
    Interpret64(I->InstrLen2,I->Instr2, registers, div32);
#if 0 //#ifndef NDEBUG
    //debug ntests values created by Interpret64 should match existing ones
    const int c = memcmp(&(I->output[i]),registers,ntests*sizeof(retval));
    if(i==0 /*c!=0*/) {
      {const int e = pthread_mutex_lock(&mutex);   assert(e==0);}
      cout<<i/npar<<" Interpret64("<<I->InstrLen2<<",I->Instr2, reg, div32);"
	  <<" ntests="<<ntests
	  <<" memcmp="<<c<<endl;
      if(c!=0) for (int j=0;j<64;j++) {
	if(j%8==0)cout<<"registers["<<j<<"]= ";
	cout<<(int)(registers[j])<<" "<<flush;
	if(j%8==7)cout<<endl;
      }
      {const int e = pthread_mutex_unlock(&mutex); assert(e==0);}
    }//endif different R0
#endif /*ndef NDEBUG*/
    memcpy(&(I->output[i]),registers,ntests*sizeof(retval));
#else /*skylake*/
    engine -> CGPengine::Interpret(I->InstrLen2,I->Instr2, i, I->output[i]);
#endif
  }
  pthread_exit(NULL);
}//end interpret

void CGPengine::CalcFitness(Individual &I)
{
#ifdef gpfunc
	//const int fd = perf_start();
#endif /*gpfunc*/
	evals++; //just for reporting stats
	I.fitness = 0;
	if(!Reg){
	  FitnessCaseNum = FitnessCase->Num();
	  const int FitnessCaseNum_ = 64*((FitnessCaseNum+63)/64);
	  //cout<<"FitnessCaseNum "<<FitnessCaseNum<<" FitnessCaseNum_ "<<FitnessCaseNum_<<endl;
	  Reg    = new retval[NumVar*FitnessCaseNum];
	  Reg64  = new retval[NumVar*FitnessCaseNum_];
	  memset(Reg64,251,NumVar*FitnessCaseNum_*sizeof(retval)); //251 not in Mackey-Glass training data
	  Output = new retval[       FitnessCaseNum];
	  for (int fcase=0;fcase<FitnessCaseNum;fcase++) {
	    for (int i=0;i<NumVar;i++) {
	      const int I = fcase*NumVar+i;
	      assert(I>=0 && I < 8*1201);
	      Reg[I] = (i<FitnessCase->Input())? FitnessCase->Input(fcase,i) : /* 0.5 */ 0;
	      const int J = fcase%npar + i*npar + (fcase/npar)*npar*NumVar;
	      assert(J>=0 && J < 8*FitnessCaseNum_);
	      assert(Reg64[J]==251);
	      Reg64[J] = Reg[I];
	    }
	    Output[fcase] = FitnessCase->Output(fcase);
	  }
#ifndef NDEBUG
	  for (int I=0;I<FitnessCaseNum_*NumVar;I++) {
	    if(I%64 ==  0) cout<<I/64/8<<" "<<(I/64)%8<<" Reg64["<<I<<"]= "<<flush;
	    cout<<(int)Reg64[I]<<" "<<flush;
	    if(I%64 == 63) cout<<endl;
	  }
#endif
	}//endif setup Reg and Reg64
#ifndef NDEBUG
	cout<<endl;
	
	for (int i=0;i<NumVar*FitnessCase->Num();i++) {
	  if(!(   ((i<128*8 && Reg[i] >=0) || Reg[i] > 27) && Reg[i] < 170))
	    cerr<<"Reg["<<i<<"]="<<(int)Reg[i]<<endl;
	  assert( ((i<128*8 && Reg[i] >=0) || Reg[i] > 27) && Reg[i] < 170);
	}
	for (int i=0;i<       FitnessCase->Num();i++)
	  assert( Output[i] > 27                        && Output[i] < 170);

	memset(I.output,253,FitnessCase->Num()*sizeof(retval)); //253 not in Mackey-Glass training data
	retval s_output[1201]; //for sanity check
	memset(s_output,252,FitnessCase->Num()*sizeof(retval)); //252 not in Mackey-Glass training data
#endif
#if (defined gpfunc) //|| (! (defined NDEBUG))
	for (int i=0;i<FitnessCaseNum;i++) {
	  Interpret(I.InstrLen, I.Instr, i, I.output[i]); //s_output[i]);
	}
	//do not break I.output[9] = 99; //break it
	//cout<<"I.output["<<9<<"]="<<(int)(I.output[9])<<endl;
	for (int i=0;i<FitnessCaseNum;i++) {
	  if(i%16==0) cout<<"I.output["<<i<<"]=";
	  cout<<(int)(I.output[i])<<" ";
	  if(i%16==15) cout<<endl;
	  //if(i>=63) break;
	}
	if(FitnessCaseNum%16 != 15) cout<<endl;
	//perf_end(fd,"CGPengine::CalcFitness(Individual) Interpret only");
#endif
#ifndef gpfunc
	/* Assume almost all overhead is here and each fitness case is similar.
	 * So give each thread a fitness case until they are all finished.
	 */
	engine = this;
	newguy = &I;
	slot = 0;
	pthread_t threads[nthreads];
	for(int i=0;i<nthreads;i++) {
	  const int e = pthread_create(&threads[i],NULL,interpret,NULL); 
	  if(e!=0){
	    cout<<"pthread_create(&threads["<<i<<"], NULL, thread_fitness) "
		<<"returned error "<<e<<endl; //At least report, perhaps its a warning try to continue
          }
        }//endfor nthreads
	/* Wait for all threads to complete */
	for(int i=0;i<nthreads;i++) {
	  const int e = pthread_join(threads[i], NULL);
	  if(e!=0){ //Report, perhaps its a warning try to continue
	    cout<<"pthread_join(threads["<<i<<"], NULL) "
		<<"returned error "<<e<<endl;
	  }
	}
#endif /*hack gpfunc*/
	for (int i=0;i<FitnessCaseNum;i++)
	{
		const int error = Output[i] - I.output[i];
		//cout<<"i"<<i<<" "<<(int)FitnessCase->Output(i)<<" "<<(int)t_output[i]<<" error="<<error<<" "
		//    <<(int)s_output[i]<<" "<<endl;
		I.fitness += error * error; //should not overflow 32bit int
#ifndef NDEBUG
#ifndef gpfunc
//		assert(I.output[i] == s_output[i]); not supported with pthreads
#endif
#endif
	}
#ifdef stats
	if(I.fitness < BestFitness) BestFitness = I.fitness;
#endif /*stats*/
}

/* replaced...
extern void GPFunc(retval v[]);
void GPFunc(retval v[]) {;};
void CGPengine::RunFitness(ofstream& out) const
{
	assert(0); //not used
	int NumVar = ::NumVar; //pickup value in GPfunc.cpp (into local copy)
	if(FitnessCase->Input() > NumVar)
		cerr<<"Warning more inputs(" <<FitnessCase->Input() <<") than registers("<<NumVar<<")\n";
	if(FitnessCase->Output() > NumVar)
		cerr<<"Warning more outputs("<<FitnessCase->Output()<<") than registers("<<NumVar<<")\n";

	double fitness = 0;
	double* f = new double[FitnessCase->Output()]; {for(int j=0; j<FitnessCase->Output();j++) f[j]=0;}
	retval* output = new retval[NumVar];
	memset(output,0,sizeof(retval)*NumVar);
	for (int i=0;i<FitnessCase->Num();i++)
	{
		{for(int j=0;j<NumVar;j++) output[j] = (j<FitnessCase->Input())? FitnessCase->Input(i,j) : /* 0.5** 0;};
		{for(int j=0;j<FitnessCase->Input();j++)  out<<FitnessCase->Input(i,j)<<" ";}
		{for(int j=0;j<FitnessCase->Output();j++) out<<FitnessCase->Output(i,j)<<" ";}

		GPFunc(output);

		{for(int j=0;j<FitnessCase->Output();j++) {
			const double error = FitnessCase->Output(i,j) - output[j];
			out<<" "<<output[j]<<" error="<<error<<" "<<flush;
			fitness += error * error;
			f[j]    += error * error;
		}}
		out<<endl;
	}
	delete[] output;

	//I.fitness += I.InstrLen; //PARSIMONY!!!

	//return Error;
	out<<"#fitness "<<fitness<<" +InstrLen ( ";
	{for(int j=0;j<FitnessCase->Output();j++) out << f[j] <<" ";}
	out<<")"<<endl;
	//GenerateCode(I);
	out<<flush;
	//out.close();

	delete[] f;
}

const double degree90 = 2*atan2(1,1);
const double pi = 2*degree90;
*/
#if 0 /*def gpfunc*/
void CGPengine::display(const int test, const retval reg[NumVar]) const //ugly or what
{
  printf("test %4d ",test);
  for(int i=0; i<NumVar;i++) printf(" %3d",reg[i]);
  printf("\n");
}
void CGPengine::display(const int test, const retval reg[NumVar], const int line, const instr &ii) const
{
  printf("test %4d ",test);
  for(int i=0; i<NumVar;i++) printf(" %3d",reg[i]);
  printf(" %3d ",line);
  fflush(stdout);
  GenerateCode(ii,cout);
  cout<<endl;
}
#else  /*gpfunc*/
void CGPengine::display(const int test, const retval reg[NumVar]) const {;}
void CGPengine::display(const int test, const retval reg[NumVar], const int line, const instr &ii) const {;}
#endif /*gpfunc*/

#define set(x) set_active(active,x)
inline void set_active(int& active, const int reg){
  assert(reg>=0 && reg < NumVar);
  active = active | (1<<reg);
  assert(active>=0 & active < (1<<NumVar));
}
#define testactive(x) test_active(active,x)
inline int test_active(int& active, const int reg){
  assert(reg>=0 && reg < NumVar);
  return (active & (1<<reg));
}
#define clear(x) clear_active(active,x)
inline void clear_active(int& active, const int reg){
  assert(reg>=0 && reg < NumVar);
  active = active & (~(1<<reg));
  assert(active>=0 & active < (1<<NumVar));
}
#define addactive(a,b,c) add_active(active,a,b,c)
inline void add_active(int& active, const int reg1, const int reg2, const int op){
  assert(reg1>=0 && reg1 < NumVar);
  //reg2 can be validvar or junk

  //check for special cases. NB protected division x/x depends on x
  if(reg1==reg2 && op==1 /*CGPengine::minus*/) return;
  set(reg1);
  if(validvar(reg2)) set(reg2);
}
   //CGPengine::Simplify(const int InstrLen, const instr *Instr, int& InstrLen2, instr *Instr2)
void CGPengine::Simplify(Individual &I, int* Needed/*=NULL*/) {
#ifdef gpfunc
  const int fd = perf_start();
#endif /*gpfunc*/
  //based active.awk r1.4
  if(!Needed) I.InstrLen2 = 0;  //as well as receiving answer Needed serves as flag to prevent changes
  int* needed = (Needed)? Needed : (int*)calloc(I.InstrLen,sizeof(int));

  int active = 0; set(0); //init_active, output in register a, ie 0
  int i;
  //backwards pass
  for(i=I.InstrLen-1; i>=0; i--){
    assert(validInstr(I.Instr[i]));
    const OP  val = I.Instr[i][0];
    const int op1 = InstrReg(I.Instr[i][1]);
    const OP  op  = I.Instr[i][2];
    const int op2 = InstrArg(I.Instr[i][3]);
    if(testactive(val)) {
      needed[i] = 1;
      clear(val);
      addactive(op1,op2,op);
    }
  }
  if(Needed) return; //dont update I
  for (i=0;i<I.InstrLen;i++) {
    if(needed[i]) {memcpy(&I.Instr2[I.InstrLen2], &I.Instr[i],sizeof(instr)); I.InstrLen2++;}
  }
  free(needed);
#ifdef gpfunc
  perf_end(fd,"Simplify");
#endif /*gpfunc*/
}
#ifdef stats
void CGPengine::Interpret(const Individual &I, retval* output2) const {
  for (int i=0;i<FitnessCaseNum;i++) {
    const size_t k = ((size_t)i)*((size_t)I.InstrLen)*((size_t)NumVar);
    assert(k>=0 && k+I.InstrLen*NumVar<=((size_t)I.InstrLen)*fcn_nv);
    assert(k>=0);
    assert(k*sizeof(retval) < output2_size);
    assert((k+I.InstrLen*NumVar)*sizeof(retval) <= output2_size);
    retval dummy;
    Interpret(I.InstrLen, I.Instr, i, dummy, &output2[k]);
  }
}
#endif /*stats*/
void CGPengine::Interpret(const int InstrLen, const instr *Instr, const int fcase, retval &output, retval* output2 /*= NULL*/) const
{
#if (defined stats) && (!(defined NDEBUG))
  if(output2) {
    if(output2_size <= 0) {
      cout<<"output2_size="<<output2_size<<endl;
      cerr<<"output2_size="<<output2_size<<endl;
      exit(99);
    }
    for(int i=0;i<InstrLen*NumVar;i++){//calloc should have ensured all are 0
      if(i*sizeof(retval) >= output2_size){
	cout<<"i="<<i<<" output2_size="<<output2_size<<endl;
	cerr<<"i="<<i<<" output2_size="<<output2_size<<endl;
	exit(99);
      }
      if(output2[i]){
	cout<<"output2["<<i<<"]="<<(int)(output2[i])<<endl;
	cerr<<"output2["<<i<<"]="<<(int)(output2[i])<<endl;
	exit(99);
      }}
  }
#endif /*stats and NDEBUG*/
#ifdef gpfunc
	const int fd = perf_start();
#endif /*gpfunc*/
	// registers
	retval reg[NumVar];
	const int I = fcase*8;
	assert(NumVar == 8);
	assert(fcase >= 0 && fcase < 1201);
	assert(I>=0 && I+7 < 8*1201);
	assert(I>=0 && I+NumVar-1 < 9728);
	memcpy(reg,&Reg[I],NumVar*sizeof(retval));

	for (int i=0;i<NumVar;i++) assert(reg[i] == FitnessCase->Input(fcase,i));
	display(fcase,reg);

	// run program
	for (int i=0;i<InstrLen;i++)
	{		
		assert(validInstr(Instr[i]));
		retval val;
		retval op1,op2;

		op1 = InstrReg(Instr[i][1],reg);
		op2 = InstrArg(Instr[i][3],reg);

		switch(int(Instr[i][2]))
		{

			case plus: 
				val = op1+op2;
				break;
			case minus: 
				val = op1-op2;
				break;
			case mul: 				
				val = op1*op2;
				break;
			case div: 
				if (op2 == 0)
					val = 0;
				else
					val = op1/op2;
				break;
#ifdef elvis
			case _acos:
				val = (op1<-1||op1>1)? 0 : acos(op1);
				break;
			case _asin:
				val = (op1<-1||op1>1)? 0 : asin(op1);
				break;
			case _atan2:
				val = (op2==0)? ((op1>=0)? degree90 : -degree90) : atan2(op1,op2);
				break;
			case _cos:
				val = cos(op1);
				break;
			case _fabs:
				val = fabs(op1);
				break;
			case _hypot:
				val = hypot(op1,op2);
				break;
			case _sign:
				val = (op1>0)? 1 : ((op1==0)? 0 : -1);
				break;
			case _sin:
				val = sin(op1);
				break;
			case _sqrt:
				val = (op1>=0)? sqrt(op1) : -sqrt(-op1);
				break;
			case _tan:
				val = ((op1-floor(op1/pi)*pi)==degree90)? 0 : tan(op1);
				break;
#endif /*elvis*/
			default:
				cerr<<"Bad opcode "<<Instr[i][2]<<endl;
				exit(99); //out <<"Bad opcode "<<I.Instr[i][2]<<endl;
			break;
		}
		// store i register
		reg[int(Instr[i][0])] = val;
#ifdef stats
		if(output2) memcpy(&output2[i*NumVar],reg,NumVar*sizeof(retval));
#endif /*stats*/
		display(fcase,reg,i,Instr[i]);
	}

	output = reg[0];

	//out<<" o="<<output[0]<<" "<<flush;
#ifdef gpfunc
	char buf[80];
	sprintf(buf,"CGPengine::Interpret(%d,Instr,%d,%d)",InstrLen,fcase,output);
	perf_end(fd,buf);
#endif /*gpfunc*/
}

 
void CGPengine::GenerateIndividual(Individual &I)
{
#ifdef elvis
	I.InstrLen = rrand(MaxInstr/2-1)+1;
	I.InstrLen2 = 0;
#else
	I.InstrLen = rrand(MaxInitInstr)+1;
#endif
	assert(FitnessCase->Num()==1201); //Mackey-Glass
	I.output = new retval[FitnessCase->Num()];
	I.Instr  = new instr[MaxInstr];
	I.Instr2 = new instr[MaxInstr];
	I.fitness = 0;

	for (int i=0;i<I.InstrLen;i++)
		GenerateInstr(I.Instr[i]);
	
	Simplify(I);
	CalcFitness(I);
}

int CGPengine::InstrArg(const OP code) const //return code not value
{
	const int x = code-IntRangeEnd-1;
	assert(code<=IntRangeEnd || (x    >= 0 &&    x < NumVar));
#ifndef elvis
	assert(code >IntRangeEnd || (code >= 0 && code < 256)); //retval == unsigned char
#endif
#ifdef NDEBUG
	return x;
#else
	if(code>IntRangeEnd) return x;       //code for Var
	else                 return INT_MAX; //code is constant, do not use
#endif /*NDEBUG*/
}

retval CGPengine::InstrArg(const OP code, const retval reg[]) const //see also VarOrVal
{
	if(code>IntRangeEnd) return reg[InstrArg(code)]; //code for Var
	else                 return code;                //code is constant
}

int CGPengine::InstrReg(const OP code) const //return code not value
{
	const int x = code-IntRangeEnd-1;
	assert(   x >= 0 &&    x < NumVar);
	return x;
}

retval CGPengine::InstrReg(const OP code, const retval reg[]) const //see also VarOrVal
{
	const int x = InstrReg(code);
#ifndef elvis
	assert(reg[x]>=0 && reg[x]<256); //retval == unsigned char
#endif
	return reg[x];
}

OP CGPengine::GenerateReg() const // Randomly choose a var
{
	return  IntRangeEnd + rrand(NumVar) + 1; //code for Var
}

OP CGPengine::GenerateArg() const // Randomly choose a var or a const
{
	if(rrand(10)<2) return GenerateReg(); //code for Var
	else            return (IntRangeEnd-IntRangeStart)*double(rand())/RAND_MAX+IntRangeStart;
}

void CGPengine::GenerateInstr(instr I)
{
	int Error=0;

	// a := b op c   // a var, b,c constant/var, op add/sub/mul/div	
	// Choose variable
	I[0] = rrand(NumVar);
	// Radomly choose a var or a const
	I[1] = GenerateReg();
	// randomly choose operator
	I[2] = rrand(NumOp);
	// Radomly choose a var or a const
	I[3] = GenerateArg();
	/*
	assert(validvar(   I[0]));
	assert(validreg(   I[1]));
	assert(validopcode(I[2]));
	cout<<"I[3] "<<I[3]<<" "<<(I[3]>=0)<<endl<<flush;
	assert((validreg  (I[3]) || validvalue(I[3])));
	*/
	assert(validInstr(I));
}


void CGPengine::GenerateCode(const instr &i, ostream& out) const
{
		assert(validInstr(i));
		out<<Var(int(i[0]));
		out<<"=";
		switch(int(i[2])) {
#ifndef elvis
		case plus :  Var(i[1],out); out<<"+"; VarOrVal(i[3],out); break;
		case minus : Var(i[1],out); out<<"-"; VarOrVal(i[3],out); break;
		case mul :   Var(i[1],out); out<<"*"; VarOrVal(i[3],out); break;
		case div :   Var(i[1],out); out<<"/"; VarOrVal(i[3],out); break;
#else /*elvis*/
		case plus :  VarOrVal(i[1],out); out<<"+"; VarOrVal(i[3],out); break;
		case minus : VarOrVal(i[1],out); out<<"-"; VarOrVal(i[3],out); break;
		case mul :   VarOrVal(i[1],out); out<<"*"; VarOrVal(i[3],out); break;
		case div :   VarOrVal(i[1],out); out<<"/"; VarOrVal(i[3],out); break;
		case _acos :
				out<<"acos("; VarOrVal(i[1],out); out<<")";  break;
		case _asin :
				out<<"asin("; VarOrVal(i[1],out); out<<")";  break;
		case _atan2 :
				out<<"atan2(";VarOrVal(i[1],out); out<<","; VarOrVal(i[3],out); out<<")";  break;
		case _cos :
				out<< "cos("; VarOrVal(i[1],out); out<<")";  break;
		case _fabs :
				out<<"fabs("; VarOrVal(i[1],out); out<<")";  break;
		case _hypot :
				out<<"hypot(";VarOrVal(i[1],out); out<<","; VarOrVal(i[3],out); out<<")";  break;
		case _sign :
				out<<"sign("; VarOrVal(i[1],out); out<<")";  break;
		case _sin :
				out<< "sin("; VarOrVal(i[1],out); out<<")"; break;
		case _sqrt :
				out<<"sqrt("; VarOrVal(i[1],out); out<<")"; break;
		case _tan :
				out<< "tan("; VarOrVal(i[1],out); out<<")";  break;
#endif
		default:
		out<< "?"; break;
		}

/*
		VarOrVal(i[1],out);

		out<<Op(int(i[2]));

		VarOrVal(i[3],out);
		out<<";";
*/
		out<<flush;
}

void CGPengine::GenerateCode(const Individual &I, const int* needed /*=NULL*/, ostream& out) const
{	
	for (int i=0;i<I.InstrLen;i++)
	{
		char buff[11];sprintf(buff,"%7d ",i);out<<buff;
		GenerateCode(I.Instr[i],out);
		if(needed && needed[i]) out<<" N";
		out<<endl;
	}
}


void CGPengine::GenerateHead(const int cnt, ostream& out) const
{
	out << "#include <math.h>\n";

	time_t t1 = time(NULL);
	out << "//" << cnt <<" "<< ctime(&t1)<<endl;

	out << "const int NumVar = "<<NumVar<<";\n";
	out << "extern void GPFunc(retval[]);\n\n";
#ifdef elvis
	out << "const double degree90 = 2*atan2(1,1);\n";
	out << "const double pi = 2*degree90;\n\n";
#endif /*elvis*/
	out << "void GPFunc(retval v[])"<<endl;
	out << "{" << endl;
#ifdef elvis
	out << "#define sqrt(x)    (((x)>=0)? sqrt(x) : -sqrt(-(x)))"<<endl;
	out << "#define sign(x)    ((x)>0)? 1 : (((x)==0)? 0 : -1)"<<endl;
	out << "#define acos(x)    (((x)<-1||(x)>1)? 0 : acos(x))"<<endl;
	out << "#define asin(x)    (((x)<-1||(x)>1)? 0 : asin(x))"<<endl;
	out << "#define atan2(x,y) (((y)==0)? (((x)>=0)? degree90 : -degree90) : atan2((x),(y)))"<<endl;
	out << "#define tan(x)     ((((x)-floor((x)/pi)*pi)==degree90)? 0 : tan(x))"<<endl;
#endif /*elvis*/
	for (int i=0;i<NumVar;i++) out << "retval "<<Var(i)<<" = v["<<i<<"];"<<endl;
}


void CGPengine::GenerateTail(ostream& out) const
{
	for (int i=0;i<NumVar;i++) out<<"v["<<i<<"] = "<<Var(i)<<";"<<endl;
	out<<"}"<<endl<<flush;
}


void CGPengine::Var(const OP i, ostream& out) const
{	
	const int op = i-IntRangeEnd-1;
	assert(op>=0 && op < NumVar);
	out<<Var(op);
}

char *CGPengine::Var(int op) const
{	
	assert(op>=0 && op < NumVar);
	return VarList[op];
}


void CGPengine::VarOrVal(const OP i, ostream& out) const //also in Interpret
{
	char buff[4];
	if (i>IntRangeEnd)
#ifdef elvis
	{		
		out<<Var(int(i-IntRangeEnd-1));
	}	
	else out<<"("<<i<<")";
#else /* make narrower and fixed width*/
	  sprintf(buff,"%-3s",Var(int(i-IntRangeEnd-1)));
	else { 
	  assert(i>=0 && i<999);
	  sprintf(buff,"%-3d",i);
	}
	out<<buff;
#endif /*elvis*/
}

/*
char *CGPengine::Int(int i)
{
	char buf[10];	
	return itoa(i,buf,10);	
}

char *CGPengine::Op(const int i) const
{
	switch(i) 
	{
		case plus : 
				return "+";
		case minus : 
				return "-";
		case mul : 
				return "*";
		case div : 
				return "/";
		case _fabs :
				return " fabs ";
		case _sign :
				return " sign ";
		case _sin :
				return " sin  ";
		case _sqrt :
				return " sqrt ";

		default:
		return "?";
	}
}
*/


