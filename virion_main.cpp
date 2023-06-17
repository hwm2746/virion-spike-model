// BEGINLICENSE
// 
// This file is part of virion-spike-model, which is distributed under the
// MIT license, as described in the LICENSE file in the top level directory
// of this project.
// 
// Author: Wonmuk Hwang
//   
// ENDLICENSE

#define SET_EXT // Do not put extern for global vars
#include "virion.h"
//#define VMD // uncomment to write vmd file

/************************************************************************/
int main(int argc, char* argv[])
{
  if (argc!=2) { 
    cout<<"Include argument in the command line: ./te input.dat" 
	<< endl; return 1;
  }
  int irun;
  ofstream fpair,fcut;
  stringstream ss;
  
  read_input(argv[1]);
  init();
  fpair.open(oprefx+"_dist.dat");
  fcut.open(oprefx+"_rcut.dat");
  write_header(fpair);  write_header(fcut);
  fpair<<"# Each data: pair_idx nn_dist(nm)"<<endl;
  fcut<<"# Number of unique min dist pairs with distance less than rcut (nm)"<<endl;
  fcut<<"# irun count"<<endl;
  for (irun=0;irun<nrun;++irun) {
    fpair<<"# irun "<<irun<<endl;
    gen_virion();
    meas_dist();
    meas_mind(fpair);
    meas_rcut(irun,fcut);
#ifdef VMD
    if (irun>5) continue; // do not write vmd files more than 5
    ss.str(""); ss.clear();
    ss<<oprefx<<"_virion_"<<irun<<".vmd";
    write_vmd(ss.str(),dpos);
    ss.str(""); ss.clear();
    ss<<oprefx<<"_nn_"<<irun<<".vmd";
    write_vmd_pair(ss.str(),mind_pair,"nn");
    ss.str(""); ss.clear();
    ss<<oprefx<<"_rcut_"<<irun<<".vmd";
    write_vmd_pair(ss.str(),rcut_pair,"rcut");    
#endif    
  }
  fpair.close();
  //write_vmd("view_virion.vmd",dpos);
  //write_vmd_pair("view_nn.vmd",mind_pair);
  //write_vmd_pair("view_rcut.vmd",rcut_pair);    
  return 0;
}
