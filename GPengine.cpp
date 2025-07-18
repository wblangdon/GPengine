// GPengine.cpp: implementation of the CGPengine class.
//
//////////////////////////////////////////////////////////////////////
//W.B.Langdon @ cs.ucl.ac.uk 23 August 2000 Elvis Hand-Eye cordination experiment
//Changes
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

//	BestFitness = FLT_MAX;

//	BestIndividual = -1;

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
void CGPengine::GenerateBestCode(const int cnt, ostream& out) const
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


void CGPengine::GenerateCode(ostream& out) const
{
	for (int i=0;i<PopSize;i++)
	{
		/*
		// Header
		cout << "double Individual(char v[])" << endl << "{" << endl;
		cout << "double a=b=c=d=0;" << endl << "a=v[0];"  << endl;
		*/
		GenerateCode(Pop[i],out);
		// footer
		//cout << "return a;" << endl;
	}
}



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
			Crossover();		
		
			// mutate loosers??
			if (randtest(PMutation))
				Mutation();		

			// calulate fittness for loosers
			Simplify(Pop[Looser1]);
			if(Changed(Looser1,Winner1) && Changed(Looser1,Winner2)) CalcFitness(Pop[Looser1]);
			Simplify(Pop[Looser2]);
			if(Changed(Looser2,Winner1) && Changed(Looser2,Winner2)) CalcFitness(Pop[Looser2]);
		}
		// no crossover, just reproduct winners
		else
			Reproduction();
		
		cnt+=2;

		if (cnt%PopSize<2)
		{	
			GenerateBestCode(cnt,out);
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
int CGPengine::Crossover()
{
	int Error=0;
	// select randomplace in winners
	int start1 = -1;
	int start2 = -1;
	int len1   = -1;
	int len2   = -1;
	ChooseXO(Pop[Winner1].InstrLen,start1,len1,
		 Pop[Winner2].InstrLen,start2,len2);
	assert(CrossLenOk(Pop[Winner1].InstrLen,len1,len2));
	assert(CrossLenOk(Pop[Winner2].InstrLen,len2,len1));

	// Create looser2 
	Cross1(Winner1,start1,len1,Winner2,start2,len2,Looser2);

	// Create looser 1	
	Cross1(Winner2,start2,len2,Winner1,start1,len1,Looser1);

/*
	if (!f)
	{
	GenerateCode(Pop[Winner1]);	
	cout << "--------------" << endl;
	GenerateCode(Pop[Winner2]);	
	cout << "start1 " << start1 << " len1 " << len1  << " start2 " << start2 << " len2 " << len2 << endl;
	GenerateCode(Pop[Looser2]);	
	cin >> Error;
	}
//	else
//		cout << "NADA"<<endl;
*/


	return Error;
}


int CGPengine::Mutation()
{
	int Error=0;

	for (int i=0;i<4;i++)
	{
		// select instruction to mutate
		int Instr1 = rrand(Pop[Looser1].InstrLen);
		int Instr2 = rrand(Pop[Looser2].InstrLen);

		MutateInstr(Instr1,Pop[Looser1]);
		MutateInstr(Instr2,Pop[Looser2]);	
	}
	return Error;
}

int CGPengine::MutateInstr(int Instr,Individual &I)
{
	int Error=0;

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

	return Error;
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

void CGPengine::CalcFitness(Individual &I)
{
	//int Error=0;
	//ofstream out;
	//if(I.InstrLen==1 && int(I.Instr[1][2])==mul)  {
	//	out.open("data.txt",ios::app);
	//out<<"\nCalcFitness of ";
	//GenerateCode(I,out);
	//}
	evals++; //just for reporting stats
	I.fitness = 0;
#ifndef NDEBUG
	memset(I.output,253,FitnessCase->Num()*sizeof(retval)); //253 not in Mackey-Glass training data
#endif
	for (int i=0;i<FitnessCase->Num();i++)
	{
		assert(i<1201 && I.output[i]==253);
#ifndef NDEBUG
		retval output; //for sanity check
		Interpret(I.InstrLen, I.Instr,  i,   output);
#endif
		Interpret(I.InstrLen2,I.Instr2, i, I.output[i]);
		assert(I.output[i] == output);
		const int error = FitnessCase->Output(i) - I.output[i];
			//out<<"i,j "<<i<<","<<j<<" "<<FitnessCase->Output(i,j)<<" "<<(int)output[j]<<" error="<<error<<" "<<endl;
			//cout<<output[j]<<" "<<flush;
		I.fitness += error * error; //should not overflow 32bit int
#ifndef NDEBUG
		const int x = i+1;
		if(x<1201) if(I.output[x]!=253) {
		  fprintf(stderr,"i=%4d I.output[%4d] %3d\n",i,x,I.output[x]);
		  fflush(NULL);
		  exit(99);
		}
#endif
	}
//	sanity(I,FitnessCase); too much even for NDEBUG
  
//	if(FitnessCase->Num()!=0) {
//		for(int j=0;j<FitnessCase->Output();j++) {I.fitness[j] = sqrt(I.fitness[j]/FitnessCase->Num());} //RMS error
//	}
	//NO parsimony I.fitness += I.InstrLen; //PARSIMONY!!!

	//return Error;
	//out<<"I.fitness "<<I.fitness<<endl;
	//GenerateCode(I);
	//out<<flush;
	//out.close();
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
#ifdef gpfunc
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
void CGPengine::Simplify(Individual &I)
{
  //based active.awk r1.4
  I.InstrLen2 = 0;
  int* needed = (int*)calloc(I.InstrLen,sizeof(int));

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
  for (i=0;i<I.InstrLen;i++) {
    if(needed[i]) {memcpy(&I.Instr2[I.InstrLen2], &I.Instr[i],sizeof(instr)); I.InstrLen2++;}
  }
  free(needed);
}
void CGPengine::Interpret(const int InstrLen, const instr *Instr, const int fcase, retval &output) const
{
	// registers
	retval reg[NumVar];
	
	for (int i=0;i<NumVar;i++) reg[i] = (i<FitnessCase->Input())? FitnessCase->Input(fcase,i) : /* 0.5 */ 0;
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

		display(fcase,reg,i,Instr[i]);
	}

	output = reg[0];

	//out<<" o="<<output[0]<<" "<<flush;
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
		case plus :  VarOrVal(i[1],out); out<<"+"; VarOrVal(i[3],out); break;
		case minus : VarOrVal(i[1],out); out<<"-"; VarOrVal(i[3],out); break;
		case mul :   VarOrVal(i[1],out); out<<"*"; VarOrVal(i[3],out); break;
		case div :   VarOrVal(i[1],out); out<<"/"; VarOrVal(i[3],out); break;
#ifdef elvis
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
*/
		out<<";";
}

void CGPengine::GenerateCode(const Individual &I, ostream& out) const
{	
	for (int i=0;i<I.InstrLen;i++)
	{
		GenerateCode(I.Instr[i],out);
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


char *CGPengine::Var(int op) const
{	
	return VarList[op];
}


void CGPengine::VarOrVal(const OP i, ostream& out) const //also in Interpret
{
	if (i>IntRangeEnd)
	{		
		out<<Var(int(i-IntRangeEnd-1));
	}	
	else out<<"("<<i<<")";
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


