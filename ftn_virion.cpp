// BEGINLICENSE
// 
// This file is part of virion-spike-model, which is distributed under the
// MIT license, as described in the LICENSE file in the top level directory
// of this project.
// 
// Author: Wonmuk Hwang
//   
// ENDLICENSE

#include "virion.h"
/************************************************************************/
void init()
{
  double rdum;
//  radius0=0.095455; // dot radius/height 0.005; //
//  height0=0.1245;
#ifdef ZHU_NAT06 /* used the way in zhu_nat06: Fig S1 caption */
  cout<<"Following the (wrong) method in Zhu et al, Nature (2006)."<<endl;
  drand0=uniform_real_distribution<double>(-1,1);
  drand1=uniform_real_distribution<double>(0.0,2*PI);
#else /* correct way */
  gaussian=normal_distribution<double>(0.0,1.0);
#endif  
  if (seed==0) { // default: use current time
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  rng.seed(seed);
  // rcut: distance cutoff, rcut0: in radian
  dotcut=rcut+2*radius0; // convert from surf distance to center-center dist
  // Done in read_input():  dotcut/=virionR; // convert to radians
  assert(fabs(dotcut)<=1.);
  dotcut=cos(dotcut); // dot product cutoff
  if (rcut<0.) dotcut=-1.; // count all pairs
  return;
}

/************************************************************************/
void gen_virion() //mt19937 &rng,normal_distribution<double>& gaussian)
{
  int j, idot,jdot,flag,niter,maxiter;
  array<double,3> r0;
  double dpcut,rdum;
#ifdef ZHU_NAT06
  double rdum1;
  double phi,theta; // phi->theta, theta->phi in regular notation
#endif
  maxiter=(ndot>1000)?(ndot/2):1000; // suitable number of max iteration
  dpcut=cos(2*radius0); // dot product cutoff
  dpos.clear();
  for (idot=0;idot<ndot;++idot) {
    flag=0; niter=0;
    while (flag==0) { // iterate until no steric clash is found
      rdum=0.;
#ifdef ZHU_NAT06
      /* used the way in zhu_nat06: Fig S1 caption*
	 phi=arccons(h), h\in [-1,1]
	 theta\in [0,2pi]
      */
      phi=acos(drand0(rng));
      theta=drand1(rng);
      rdum1=sin(phi);
      r0[0]=rdum1*cos(theta); r0[1]=rdum1*sin(theta); r0[2]=cos(phi);
      for (j=0;j<3;j++) rdum+=r0[j]*r0[j];
#else /* use the correct way */
      for (j=0;j<3;j++) {
	r0[j]=gaussian(rng); rdum+=r0[j]*r0[j];
      }
#endif
      rdum=sqrt(rdum);
      for (j=0;j<3;j++) { r0[j]/=rdum; } // normalize
      flag=1;
      if (idot>0) { // check distance with previous dots
	for (jdot=0;jdot<idot;++jdot) {
	  rdum=dotprod(r0,dpos[jdot],3);
	  if (rdum>dpcut) { flag=0; break; } // close one found. Try again
	}
      }
      ++niter;
      if (niter>maxiter) {
	cout<<"No place to put a new dot after maximum number of iterations."
	    <<endl;
	return;
      }
    } // while (flag==0) { // iterate until no steric clash is found
    assert(flag==1);
    dpos.push_back(r0);
  } //  for (idot=0;idot<ndot;++idot) {
  return;
}

/************************************************************************/
void meas_dist()
{ // measure pairwise center-to-center distances (dot products)
  int idot,jdot;
  double rdum;
  multimap<double,int> dpmap;
  dp0.clear();

  // ndot^2 operations, but safer and simpler in meas_mind()
  for (idot=0;idot<ndot;++idot) {  
    dpmap.clear();
    for (jdot=0;jdot<ndot;++jdot) {
      if (idot==jdot) continue;
      rdum=dotprod(dpos[idot],dpos[jdot],3);
      dpmap.insert(pair<double,int>(rdum,jdot));
    } // for (jdot=idot;jdot<ndot;++jdot) {

//  for (idot=0;idot<(ndot-1);++idot) { // ndot*(ndot-1)/2 operations
//    dpmap.clear();
//    for (jdot=idot+1;jdot<ndot;++jdot) {
//      rdum=dotprod(dpos[idot],dpos[jdot],3);
//      dpmap.insert(pair<double,int>(rdum,jdot));
//    } // for (jdot=idot;jdot<ndot;++jdot) {
    dp0.push_back(dpmap);
  } // for (idot=0;idot<(ndot-1);++idot) {
}
  
/************************************************************************/
void meas_mind(ofstream &fpair)
{ /* find minimum surf-to-surf distance pairs. Since dp0[idot] is a
    multimap, the first element is the closest.
   */
  int idot,jdot,jmin,flag,npair;
  double rdum,dpmax; // max dot product (minimum distance)
  //  multimap<int,int> mind_pair; // moved to virion.h
  multimap<int,int>::iterator it;
  multimap<double,int>::reverse_iterator rjt;
  
  mind_pair.clear();
  npair=0;
  //fpair<<"# irun= "<<irun<<endl;
  for (idot=0;idot<ndot;++idot) {
    dpmax=dp0[idot].rbegin()->first; // min dist for jdot>idot
    jmin=dp0[idot].rbegin()->second;
    if (jmin<idot) {
      it=mind_pair.find(jmin);
      if (it!=mind_pair.end()) {
	if (it->second ==idot) continue; // (jmin,idot) already registered
      }
    }
    mind_pair.insert(pair<int,int>(idot,jmin));
    // surface-to-surface distance between dots
    fpair<<setw(3)<<npair<<" "<<setw(6)<<setprecision(4)<<fixed
	 <<virionR*(acos(dpmax)-2.*radius0)<<"  "<<idot<<" "<<jmin<<endl;
    ++npair;
  }
  return;
}

/************************************************************************/
void meas_rcut(int irun, ofstream &fcut)
{ /* Find and count the number of unique pairs that are within distance rcut.
   */
  int idot,jdot,jmin,flag,npair;
  double rdum;
  multimap<int,int>::iterator it1,it2;
  pair<multimap<int,int>::iterator,multimap<int,int>::iterator> itrng;
  multimap<double,int>::reverse_iterator rit;
  
  rcut_pair.clear();
  npair=0;
  for (idot=0;idot<ndot;++idot) {
    for (rit=dp0[idot].rbegin();rit!=dp0[idot].rend();++rit) {
      rdum=rit->first;
      if (rdum>dotcut) { // closer distance -> greater dot product
	jdot=rit->second;
	if (jdot<idot) {
	  it1=rcut_pair.find(jdot);
	  if (it1!=rcut_pair.end()) {
	    itrng=rcut_pair.equal_range(jdot);
	    flag=0;
	    for (it2=itrng.first;it2!=itrng.second;++it2) {
	      if (it2->second == idot) { // (jdot,idot) already registered
		flag=1; break;
	      }
	    }
	    if (flag==1) continue; 
	  }
	}
	rcut_pair.insert(pair<int,int>(idot,jdot));
	++npair;
      } // if (rdum>rcut) {
      else break;
    } // for (rit=dp0[idot].rbegin();rit!=dp0[idot].rend();+rit) {
  }
  fcut<<setw(6)<<irun<<" "<<setw(6)<<npair<<endl;

  //ofstream ffo;
  //// For debug 0) (DO NOT DELETE)
  //ffo.open("temp0.dat"); // for debug
  //jmin=0; // use jmin as a counter
  //for (it=mind_pair.begin();it!=mind_pair.end();++it) {
  //  ffo<<setw(2)<<it->first<<" "<<it->second<<"  ";
  //  if (++jmin %5 ==0) ffo<<endl;
  //}
  //ffo<<endl;
  //ffo.close();

  return;
}

/************************************************************************/
void read_input(string ifname)
{
  ifstream fin(ifname);
  stringstream ss;
  string sdum,sdum1;
  // default values
  nm_R=5.25; nm_H=13.7; virionR=55;
  //radius0=0.095455;   height0=0.1245; // default spike dimensions (for hiv)
  seed=0; ndot=14; nrun=1; nm_rcut=-1.; oprefx="temp"; // default
  while (!fin.eof()) {
    getline(fin,sdum);
    sdum1.clear();
    sdum1=sdum.substr(0,sdum.find("#")); // remove comment
    ss.str(""); ss.clear(); ss<< sdum1;
    while (ss.good()) {
      sdum.clear();
      ss >> sdum;
      if (sdum=="") continue; // empty command
      else if (sdum=="ndot") ss>>ndot;
      else if (sdum=="nrun") ss>>nrun;
      else if (sdum=="oprefx") ss>>oprefx;
      else if (sdum=="rcut") ss>>nm_rcut;
      else if (sdum=="spike_radius") ss>>nm_R;
      else if (sdum=="spike_height") ss>>nm_H;
      else if (sdum=="seed") ss>>seed; // should be <=10 digits
      
      else if (sdum=="virionR") ss>>virionR;
      else {
	cout<<"Unrecognized token: "<<sdum<<endl;
      }
      break; // process the next command
    }
  }
  // scale lengths
  radius0=nm_R/virionR; height0=nm_H/virionR;
  rcut=(nm_rcut<0.)?-1.:(nm_rcut/virionR); // in radians
  return;
}

/************************************************************************/
void write_header(ofstream &fout)
{
  fout<<"# ndot= "<<ndot<<endl;
  fout<<"# nrun= "<<nrun<<endl;
  fout<<"# seed= "<<seed<<endl;
  fout<<"# spike radius/height (nm)= "<<nm_R<<" "<<nm_H<<endl;
  fout<<"# rcut (nm)= "<<nm_rcut<<endl;
  fout<<"# virionR (nm)= "<<virionR<<endl;
  return;
}
/************************************************************************/
void write_vmd(string ofname,vector<array<double,3>> &dpos)
{
  int j,idot;
  // make tip height after adding hemisphere as height0
  double height1=1.+height0-0.5*radius0; 
  array<double,3> r1;
  ofstream ffo(ofname);
  ffo<<"# Output of virion_main.cpp"<<endl;
  ffo<<"# rng seed: "<<seed<<endl;
  ffo<<"# Number of dots: "<<ndot<<"  radius: "<<radius0
     <<"  height: "<<height0<<endl;
  ffo<<"mol new"<<endl;
  ffo<<"graphics top color blue"<<endl;
  ffo<<"graphics top material Transparent"<<endl;
  ffo<<"graphics top sphere {0 0 0} radius 1 resolution 256"<<endl;
  ffo<<"mol rename top sphere"<<endl;

  ffo<<"mol new"<<endl;
  ffo<<"graphics top color red"<<endl;
  for (idot=0;idot<ndot;++idot) {
    for (j=0;j<3;j++) { r1[j]=dpos[idot][j]*height1; } 
    ffo<<"graphics top sphere { ";
    for (j=0;j<3;j++) {
      ffo<<setw(10)<<setprecision(6)<<r1[j]<<" ";
    }
    ffo<<"} radius "<<radius0<<" resolution 64"<<endl;
    ffo<<"graphics top cylinder { ";
    for (j=0;j<3;j++) {
      ffo<<setw(10)<<setprecision(6)<<dpos[idot][j]<<" ";
    }
    ffo<<" } { ";
    for (j=0;j<3;j++) {
      ffo<<setw(10)<<setprecision(6)<<r1[j]<<" ";
    }
    ffo<<"} radius "<<radius0<<" resolution 64 filled yes"<<endl;
  } // for (idot=0;idot<ndot;++idot) {
  ffo<<"mol rename top spikes"<<endl;
  ffo<<"color Display Background white"<<endl;
  ffo.close();
  return;
}

/************************************************************************/
void write_vmd_pair(string ofname,multimap<int,int> &mpair,string molname)
{ // Generate vmd file showing mpair
  ofstream ffo(ofname);
  int j,idot,jdot;
  double height1=1.+height0-0.5*radius0;
  array<double,3> r0,r1;
  multimap<int,int>::iterator it;

  ffo<<"# output of ftn_virion.cpp: Check nearest neighbor"<<endl;
  ffo<<"mol new"<<endl;
  ffo<<"graphics top color yellow"<<endl;
  //k=0;
  for (it=mpair.begin();it!=mpair.end();++it) {
    //ffo<<"graphics top color colorID "<<k<<endl;
    idot=it->first; jdot=it->second;
    for (j=0;j<3;j++) { r0[j]=dpos[idot][j]*height1; } 
    for (j=0;j<3;j++) { r1[j]=dpos[jdot][j]*height1; } 
    ffo<<"graphics top cylinder { ";
    for (j=0;j<3;j++) {
      ffo<<setw(10)<<setprecision(6)<<r0[j]<<" ";
    }
    ffo<<" } { ";
    for (j=0;j<3;j++) {
      ffo<<setw(10)<<setprecision(6)<<r1[j]<<" ";
    }
    ffo<<"} radius "<<radius0*0.5<<" resolution 64 filled yes"<<endl;
    //++k;
  }
  ffo<<"mol rename top "<<molname<<endl;
  ffo.close();
  return;
}
