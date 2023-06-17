// BEGINLICENSE
// 
// This file is part of virion-spike-model, which is distributed under the
// MIT license, as described in the LICENSE file in the top level directory
// of this project.
// 
// Author: Wonmuk Hwang
//   
// ENDLICENSE

#ifndef VIRION
#define VIRION

#ifdef SET_EXT
  #define EXT_NAME 
#else
  #define EXT_NAME extern
#endif

//#define ZHU_NAT06 // uncomment to use Zhu et al, Nature (2006) SI Fig 6

#include <iostream>
#include <iomanip>      // std::setprecision
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <cassert>
#include <random>
#include <chrono>
using namespace std;
#define PI 3.14159265358979

EXT_NAME int ndot,nrun; // ndot: number of dots (gp160) 1e5; // 73: SIV, 14: HIV
EXT_NAME double radius0,height0; // dot radius/height scaled by virionR
EXT_NAME double nm_R,nm_H; // dot radius/height in nm
// rcut,dotcut: for counting unique pairs that are within rcut
EXT_NAME double virionR, rcut,nm_rcut,dotcut; // virionR: virion radius
// nm_rcut: rcut in nm unit, rcut: scaled by virionR
EXT_NAME unsigned seed;
EXT_NAME vector<array<double,3>> dpos; // dot position
// dp0[idot]: multimap of <dot product,jdot> with other dots
EXT_NAME vector<multimap<double,int>> dp0;
EXT_NAME multimap<int,int> mind_pair,rcut_pair; // nearest neighbor pair indices
EXT_NAME mt19937 rng;
#ifdef ZHU_NAT06 /* used the way in zhu_nat06: Fig S1 caption */
EXT_NAME uniform_real_distribution<double> drand0,drand1; 
#else
EXT_NAME normal_distribution<double> gaussian; 
#endif
EXT_NAME string oprefx;

/************************************************************************/
// vector_stl.cpp
double dotprod(array<double,3>& a, array<double,3>& b, int N);
void crossprod(array<double,3>& a, array<double,3>& b, array<double,3>& c);
double getAngle(array<double,3> &bvec1, array<double,3> &bvec2, char flag);
void matmul3d(array<array<double,3>,3>& a, array<array<double,3>,3>& b,
	      array<array<double,3>,3>& c);
double normalizeVec(array<double,3> &a, int N);
void rotate_dir3d(array<double,3>& a, array<double,3>& u, array<double,3>& c,
		  double theta);
int levi_civita(int i, int j, int k);
void rotmat_dir3d(array<double,3>& u, double theta,
		  array<array<double,3>,3>& r0);

/************************************************************************/
// ftn_virion.cpp
void init();
void gen_virion(); //mt19937 &rng,normal_distribution<double>& gaussian);
void meas_dist();
void meas_mind(ofstream &fpair);
void meas_rcut(int irun, ofstream &fcut);
void read_input(string ifname);
void write_header(ofstream &fout);
void write_vmd(string ofname,vector<array<double,3>> &dpos);
void write_vmd_pair(string ofname,multimap<int,int> &mpair,string molname);

#endif

