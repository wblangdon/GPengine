//main.cpp of GPengine
//W.B.Langdon @ cs.ucl.ac.uk 23 August 2000 Elvis Hand-Eye cordination experiment
//Changes
//WBL 30 Dec 2025 display GenerateLimit
//WBL 21 Jul 2025 Make log optional, add pthreads
//WBL 19 Jun 2025 Add seed and display last pseudo-random number
//WBL 19 Jun 2025 Redo Mackey-Glass cf geccolb.tex r1.20, complex-systems.tex r1.30
//WBL 19 Jun 2025 For Linux g++ 11.5.0
//WBL 23 Aug 00 Add input file specification class, change FitnessCase
//WBL 23 Aug 00 TODO Add command line processing

/* compile:
  g++ -o GPrun GPengine.o -fpermissive -fmax-errors=2 main.cpp 
 */

//#include <stdlib.h>
//#include <iostream.h>

#include "GPengine.h"
//#include <assert.h>
#include <sys/timeb.h>
#include <sys/sysinfo.h>

BOOL dataspec::add(const int pos) 
{
	//printf("add %d\n",pos);
	if(pos<=0) return false;
	chain* prev = first;
	chain* p = first; for(; p!=NULL;) {
		if(p->pos==pos) return false;
		prev = p;
		p = p->next;
	}
	p = new chain;
	p->pos  = pos;
	p->next = NULL;
	if(first==NULL) first = p;
	else       prev->next = p;
	size++;
	if(pos > max) max = pos;
	return true;
}

dataspec::~dataspec() 
{
	for(chain* p = first; p!=NULL;) {
		chain* prev = p;
		p = p->next;
		delete prev;
	}
	first = NULL;
}

int dataspec::Pos(const int i) const {
	assert(0<=i && i<size && first!=NULL);
	int j = i;
	chain* p = first; for(; j>0;j--) {
		p = p->next;
		assert(p!=NULL);
	}
	return p->pos;
}

int dataspec::Size() const { return size; }	
int dataspec::Max()  const { return max; }

//int writefunction(char *filename);

#ifdef gpfunc
#include "GPfunc.cpp"
#endif

int nthreads = -1;
extern const char* rev;
extern const char* rev2;
//Bit of a klundge cannot set nthreads to zero, eg T0 fails
int main(const int argc, const char *argv[])
{
    if (argc < 5)    /* Test for correct number of arguments */
    {
        fprintf(stderr, "Usage: %s <infile> [outfile] <innn>* <onnn>* [sseed] [tthreads]\n", argv[0]);
        exit(1);
    }

	dataspec input;
	dataspec output;
	int seed = -1;
	const int ofile = (strchr(argv[2],'.')!=NULL); //eg file.out
	for(int i = ofile? 3:2; i<argc;i++) {
		if(((*argv[i]=='i' || *argv[i]=='I') && ( input.add(atoi(&argv[i][1]))))  ||
		   ((*argv[i]=='o' || *argv[i]=='O') && (output.add(atoi(&argv[i][1]))))  ||
		   ((*argv[i]=='s' || *argv[i]=='S') &&       (seed=atoi(&argv[i][1])))   ||
		   ((*argv[i]=='t' || *argv[i]=='T') &&   (nthreads=atoi(&argv[i][1]))) ) {}
		else {
			fprintf(stderr, "Bad Column specifier %s. Usage: innn or onnn \n", argv[i]);
	        exit(1);
		}
	}

	ofstream log;
	if(ofile) log.open(argv[2], ios::app);
#ifdef elvis
	cout<<  "#GPengine $Revision 1.00 $ WBL 29 August 2000 ";
	if(log) log << "//GPengine $Revision 1.00 $ WBL 29 August 2000 ";
#else
	cout<<  "#GPengine $Revision: 1.25 $";
	if(rev  && strlen(rev) >11) cout<<" rev="<<&rev[11]<<flush;
	if(rev2 && strlen(rev2)>11) cout<<" re2="<<&rev2[11]<<flush;
	cout<<  " WBL December 2025 "<<flush;
#endif
	{for(int i=0; i<argc; i++) cout<< argv[i] << " ";}
	{for(int i=0; i<argc; i++) if(log) log << argv[i] << " ";}
	if(seed>0) {
	  cout<< "seed=" <<seed<< " ";
	  if(log) log << "seed=" <<seed<< " ";
	  srand(seed);
	}
        { const int n = get_nprocs();
	  if(nthreads == -1) nthreads = n;
	  else if(nthreads > n)
	    cerr << "Warning asked for more threads "
		 << nthreads << " than we have "<<n<<endl;
	}
	cout<< "threads=" <<nthreads<< " ";
	if(log) log << "threads=" <<nthreads<< " ";
	if(nthreads < 1){
	  cerr << "Bad nthreads "<<nthreads<<endl;
	  exit(1);
	}
	if(log) {
	  cout<< argv[2] << " ";
	  log << argv[2] << " ";
	}
	cout<< "GenerateLimit=" <<GenerateLimit<< " ";
	if(log) cout<< "GenerateLimit=" <<GenerateLimit<< " ";
	//https://stackoverflow.com/questions/36927699/how-to-correctly-use-ctime-to-print-different-time-stamps
	{const time_t t1 = time(NULL);
	cout<< ctime(&t1) << endl;
	if(log) log << ctime(&t1) << endl;
	}
	CGPengine GP;
	
	GP.Init();

	GP.LoadTrainingData(argv[1],input,output);

	//GP.write(argv[2]);

#ifdef gpfunc
	GP.LoadRun();
#else
	GP.GeneratePop();
	GP.GenerateBestCode( 0/*,log*/);

	GP.Evolve(log);

	GP.GenerateBestCode(-1/*,log*/);
#ifndef stats
	if(log) GP.GenerateBestCode(log);
#endif /*not stats*/

	const int last_rand = rand();
	cout<<  "#last random number "<<last_rand<<" ";
	if(log) log << "//last random number "<<last_rand<<" ";
	{const time_t t1 = time(NULL);
	cout<< ctime(&t1) << endl;
	if(log) log << ctime(&t1) << endl;
	}
	if(log) log.close();
#endif /*gpfunc*/

	return 0;
}

/*
int f(int x)
{

	return (x*123*x);

	int sum=1;
	for (int i=1;i<=x;i++)
		sum +=i;

	return sum;

}

**
void _FitnessCase::write(const char *filename) const
{
	FILE *file;

	file = fopen(filename, "w+" );
	
	fprintf(file,"#%d\n",Num());
	for (int i=0;i<Num();i++) {
		{for(int j=0; j<Input();  j++) fprintf( file,"%f ",(float)Input(i,j)); }		
		{for(int j=0; j<Output(); j++) fprintf( file,"%f ",(float)Output(i,j));}
		fprintf(file,"\n");
	}
	fclose(file);
}
*/

