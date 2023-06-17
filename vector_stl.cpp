// BEGINLICENSE
// 
// This file is part of virion-spike-model, which is distributed under the
// MIT license, as described in the LICENSE file in the top level directory
// of this project.
// 
// Author: Wonmuk Hwang
//   
// ENDLICENSE

/* vector_stl.cpp: stl version of functions in vectors.cpp. 
  Compile with -fconcepts
*/

#include <iostream>
#include <iomanip>      // std::setprecision
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
using namespace std ;

#define TOL 1e-8
#define PI 3.14159265358979

/************************************************************************/
double dotprod(array<double,3>& a, array<double,3>& b, int N)
{ // Calculates \vec a\cdot \vec b in 3d. N is redundant.
  // Note that this routine does not check if a & b are of the same size
  double r=0.0;  
  for (int i=0;i<N;i++) r+=a[i]*b[i];
  return r;
}

/************************************************************************/
void crossprod(array<double,3>& a, array<double,3>& b, array<double,3>& c)
{ // Calculates \vec c= \vec a\times \vec b in 3-dimension
  c[0]=a[1]*b[2]-a[2]*b[1]; 
  c[1]=a[2]*b[0]-a[0]*b[2]; 
  c[2]=a[0]*b[1]-a[1]*b[0]; 
} 

/************************************************************************/
double getAngle(array<double,3> &bvec1, array<double,3> &bvec2, char flag)
{ /* Find angle between bvec1 & bvec2. 
     flag={0,1}: 2-Dim angle. 
        flag=0: Range is (0,PI)
        flag=1: Range is (-PI,PI). The sign is defined in the
           counterclockwise manner from bvec1 to bvec2.  
     flag=3: 3-dim angle
 */
  int ndim;
  double r1, r2, r3, phi;
  ndim=(flag<3)?2:3;
  r1=dotprod(bvec1,bvec1,ndim);  r2=dotprod(bvec2,bvec2,ndim);
  r1=sqrt(r1); r2=sqrt(r2);
  r3=dotprod(bvec1,bvec2,ndim)/(r1*r2);
  r3=(r3>1.0)?(1.0-TOL):r3; // take care of round off error
  r3=(r3<-1.0)?(-1.0+TOL):r3;
  phi=acos(r3);
  if (flag==1) { // get the sign 
    r3=bvec1[0]*bvec2[1]-bvec1[1]*bvec2[0]; // cross product
    if (r3<0) phi*=-1.0;
  }
  return phi;
}

/************************************************************************/
void matmul3d(array<array<double,3>,3>& a, array<array<double,3>,3>& b,
	      array<array<double,3>,3>& c)
{ // carry out axb=c
  int i,j,k;
  for (i=0;i<3;i++) { for (j=0;j<3;j++) c[i][j]=0.; }
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) c[i][j]+=a[i][k]*b[k][j];
    }
  }
  return;
}

///************************************************************************/
//double normalizeVec(auto &a, int N)
//{ // normalize vector a in N-dim. returns norm of a[] before normalization
//  double rdum=sqrt(dotprod(a,a,N));
//  for (auto& r : a ) r/=rdum;
//  //  for (int i=0;i<N;i++) a[i]/=rdum;
//  return rdum;
//}

/************************************************************************/
double normalizeVec(array<double,3> &a, int N)
{ // normalize vector a in 3D. returns norm of a[] before normalization
  double rdum=sqrt(dotprod(a,a,N)); // N is redundant.
  for (int i=0;i<N;i++) a[i]/=rdum;
  return rdum;
}

/************************************************************************/
void rotate_dir3d(array<double,3>& a, array<double,3>& u, array<double,3>& c,
		  double theta)
{ /* Rotate vector a about direction specified by unit vector u by angle
     theta in 3D. See doc/rotate_vector210515.pdf for explanation.
     Store rotated vector to c[].

     NOTE: Put this function below dotprod() and crossprod() in the code.
  */
  int i;
  double rdum,c0,s0;
  array<double,3> p;
  // check normalization of u
  rdum=dotprod(u,u,3);
  if (fabs(rdum-1)>1.0e-5) {
    cout<<" Direction vector not normalized. Normalizing here."<<endl;
    rdum=sqrt(rdum);
    for (i=0;i<3;i++) u[i]/=rdum;
  }
  rdum=dotprod(u,a,3);
  crossprod(u,a,p); c0=cos(theta); s0=sin(theta);
  rdum=(1.-c0)*rdum;
  for (i=0; i<3;i++) {
    c[i]=c0*a[i]+rdum*u[i]+s0*p[i];
  }
  return;
}

/************************************************************************/
array<double,3>  subtract(array<double,3>& a, array<double,3>& b, int N)
{ // return a-b
  array<double,3> adum;
  for (int i=0;i<3;i++) adum[i]=a[i]-b[i];
  return adum;
}

/************************************************************************/
int levi_civita(int i, int j, int k)
{ // Levi-Civita symbol. i,j,k=0,1,2
  if ((i<0)||(i>2)) {
    cout<<"ERROR: invalid first argument: "<<i<<endl;
    exit(-1);
  }
  if ((j<0)||(j>2)) {
    cout<<"ERROR: invalid second argument: "<<j<<endl;
    exit(-1);
  }
  if ((k<0)||(k>2)) {
    cout<<"ERROR: invalid second argument: "<<k<<endl;
    exit(-1);
  }
  if ((i==j)||(i==k)||(j==k)) return 0;
  else {
    if ((i==0)&&(j==1)&&(k==2)) return 1;
    else if ((i==1)&&(j==2)&&(k==0)) return 1;
    else if ((i==2)&&(j==0)&&(k==1)) return 1;
    else return -1;
  }
}

/************************************************************************/
void rotmat_dir3d(array<double,3>& u, double theta,
		  array<array<double,3>,3>& r0)
{ /* Given axis of roation u and rotation angle theta, build 3x3 rotation
     matrix r0. r0 should rotate a vector in the same way as rotarte_dir3d()
     does. r0_{ij} is:

    u_i*u_j*(1-cos(theta))+\delta_{ij}*cos(theta)+\epsilon_{ikj}*u_k*sin(theta)

    where \delta_{ij}: Kronecker-delta,
          \epsilon_{ikj}: Levi-Civita symbol.

    See my notes on quaternion for details.
*/
  int i,j,k;
  double rdum,c0,s0;
  array<double,3> p;
  // check normalization of u
  rdum=dotprod(u,u,3);
  if (fabs(rdum-1)>1.0e-5) {
    cout<<" Axis of rotation not normalized. Normalizing here."<<endl;
    rdum=sqrt(rdum);
    for (i=0;i<3;i++) u[i]/=rdum;
  }
  c0=cos(theta); s0=sin(theta);
  rdum=(1.-c0);
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) r0[i][j]=rdum*u[i]*u[j]; // u_i*u_j*(1-cos(theta))
    r0[i][i]+=c0; // +\delta_{ij}*cos(theta)
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) r0[i][j]+=(double)levi_civita(i,k,j)*u[k]*s0;
    }
  }
  return;
}

