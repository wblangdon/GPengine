//main.cpp of GPengine
//W.B.Langdon @ cs.ucl.ac.uk 23 August 2000 Elvis Hand-Eye cordination experiment
//Changes
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

int main(const int argc, const char *argv[])
{
    if (argc < 5)    /* Test for correct number of arguments */
    {
        fprintf(stderr, "Usage: %s <infile> <outfile> <innn>* <onnn>* seed\n", argv[0]);
        exit(1);
    }

	dataspec input;
	dataspec output;
	int seed = -1;
	for(int i=3;i<argc;i++) {
		if(((*argv[i]=='i' || *argv[i]=='I') && ( input.add(atoi(&argv[i][1]))))  ||
		   ((*argv[i]=='o' || *argv[i]=='O') && (output.add(atoi(&argv[i][1]))))  ||
		   ((*argv[i]=='s' || *argv[i]=='S') &&       (seed=atoi(&argv[i][1]))) ) {}
		else {
			fprintf(stderr, "Bad Column specifier %s. Usage: innn or onnn \n", argv[i]);
	        exit(1);
		}
	}

	ofstream log;
	log.open(argv[2], ios::app);
#ifdef elvis
	cout<<  "#GPengine $Revision 1.00 $ WBL 29 August 2000 ";
	log << "//GPengine $Revision 1.00 $ WBL 29 August 2000 ";
#else
	cout<<  "#GPengine $Revision: 1.17 $ WBL June 2025 ";
	log << "//GPengine $Revision: 1.17 $ WBL June 2025 ";
#endif
	{for(int i=0; i<argc; i++) cout<< argv[i] << " ";}
	{for(int i=0; i<argc; i++) log << argv[i] << " ";}
	if(seed>0) {
	  cout<< "seed " <<seed<< " ";
	  log << "seed " <<seed<< " ";
	  srand(seed);
	}
	//https://stackoverflow.com/questions/36927699/how-to-correctly-use-ctime-to-print-different-time-stamps
	time_t t1 = time(NULL);
	cout<< ctime(&t1) << endl;
	log << ctime(&t1) << endl;

	CGPengine GP;	
	
	GP.Init();

	GP.LoadTrainingData(argv[1],input,output);

	//GP.write(argv[2]);

#ifdef gpfunc
	GP.LoadRun();
#else
	GP.GeneratePop();
	GP.GenerateBestCode( 0,log);

	GP.Evolve(log);

	GP.GenerateBestCode(-1,log);

	//getch();

	const int last_rand = rand();
	cout<<  "#last random number "<<last_rand<<endl;
	log << "//last random number "<<last_rand<<endl;
	log.close();
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

