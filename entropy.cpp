/*WBL 11 August 2025*/
#ifndef NDEBUG
extern const char* rev2 = "$Revision: 1.11 $";
#else
extern const char* rev2 = "$Revision: 1.11 $ NDEBUG";
#endif

//Modifications:
//WBL 30 Dec 2025 make cout depend on NDEBUG

/*
 * using STL for hash map so make this a separate file
 *rather than part of GPengine.cpp
 *
 *Mackey-Glass has 8 inputs each of a byte
  so all will fit into unsigned 64 bit int
*/

/* Increment 
https://stackoverflow.com/questions/56390374/how-to-increment-a-particular-value-of-a-hash-table-without-changing-its-key
*/

//Extend Microsoft.com example
//https://learn.microsoft.com/en-us/cpp/standard-library/unordered-map-class?view=msvc-170#example

#include <cmath>
#include <unordered_map>
#include <iostream>
using namespace std;

typedef unordered_map<uint64_t, int> Mymap;
Mymap dat;
int dat_n = 0;

void clear(){
  dat_n = 0;
  dat.clear();
}

void increment(const uint64_t key){
  dat_n    += 1;
  dat[key] += 1;
}

void hex8(const uint64_t key,ostream& out) {
  uint64_t t = key;
  for(int i=0;i<8;i++){
    char buff[80];
    sprintf(buff,"%d,",t & 255);
    out<<buff<<flush;
    t = t>>8;
  }
}

inline double log2(double x){ return log(x)/log(2.0);}
double entropy(){
  if(dat_n != 1201) {
    cerr<<"entropy() bad dat_n="<<dat_n<<endl;
    exit(99);
  }
  if(dat_n == 0) return 0.0;
  double h = 0;
  for (Mymap::const_iterator it = dat.begin();
       it != dat.end(); ++it) {
    //cout << " [ " << it->first << " " << it->second << " ]";
#ifndef NDEBUG
    hex8(it->first,cout); cout << ":" << it->second << " ";
#endif
    const double p = double(it->second)/double(dat_n);
    const double h_ = -p * log2(p);
    //cout << " p == " << p << " h_ == " << h_ <<endl;
    h += h_;
  }
#ifndef NDEBUG
  cout << endl;
#endif
  return h;
}

/*
int main()
{
    const uint64_t a  = 1;
    const uint64_t b  = 0x736f6d6570736575;
    const uint64_t c  = 0xf36f6d6570736575;
    const uint64_t d  = 0;
    //increment(a);
    increment(b);
    increment(b);
    for(int i=0;i<2;i++) increment(c);
  			
    // find and show elements
    const double h1 = entropy();
    cout <<"entropy() " << h1 << endl;
    clear();
    for(int i=0;i<20;i++) increment(a);
    for(int i=0;i<20;i++) increment(b);
    for(int i=0;i<20;i++) increment(c);
    for(int i=0;i<20;i++) increment(d);
    const double h = entropy();
    cout <<"entropy() == " << h << endl;
    
    return (0);
}
*/
