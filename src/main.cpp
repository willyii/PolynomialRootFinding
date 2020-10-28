// #include <math.h>

// #include <iostream>
// #include <vector>

#include "test.h"

#include "isolate.h"

using namespace std;



int main(int argc, char const *argv[]) {
  /* code */
  // testBudan();
  

  vector<double> tmp = {1.0, -1.0};
  Poly test1 = Poly(tmp);
  test1 = test1.pow(2);
  cout<<test1<<endl;
  // Budan util;
  // test1 = util.timeToP(test1, -2);
  // cout<<test1<<endl;
  // // tmp = {0.0, 4.0};
  // // Poly test2 = Poly(tmp);
  // Poly test = test1;// * test2 * test2;


  // Isolate tttt = Isolate();
  // vector<Poly> ans = tttt.squareFreeDecompo(test);
  // cout<<"Square Free Decompo ans: "<<endl;
  // for(auto t : ans )
  // cout<<t<<"\n";

  return 0;
}
