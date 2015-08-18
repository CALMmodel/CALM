/* relativistic algebra: implementation */

#include "specrel.h"
#include <math.h>
#include <iostream>
using namespace std;

vector4::vector4()
{
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;
  v[3] = 0.;
}

vector4::vector4(double a0, double a1, double a2, double a3)
{
  v[0] = a0;
  v[1] = a1;
  v[2] = a2;
  v[3] = a3;
}

vector4::vector4(double a[4])
{
  v[0] = a[0];
  v[1] = a[1];
  v[2] = a[2];
  v[3] = a[3];
}

vector4::vector4(const vector4 & rhs)
{
  v[0] = rhs.v[0];
  v[1] = rhs.v[1];
  v[2] = rhs.v[2];
  v[3] = rhs.v[3];
}

vector4::~vector4()
{ }

double  &vector4::operator[](int index)
{
  if ((index<0)||(index>=4)){
    cout << "\n\a ****  4-vector index under/overflow  ****\n\n";
    return v[0];
  }
  return v[index];
}

vector4 operator+(const vector4 & v1, const vector4 & v2)
{
  vector4 hv;
  hv[0] = v1.v[0] + v2.v[0];
  hv[1] = v1.v[1] + v2.v[1];
  hv[2] = v1.v[2] + v2.v[2];
  hv[3] = v1.v[3] + v2.v[3];
  return hv;
}

vector4 operator-(const vector4 & v1, const vector4 & v2)
{
  vector4 hv;
  hv[0] = v1.v[0] - v2.v[0];
  hv[1] = v1.v[1] - v2.v[1];
  hv[2] = v1.v[2] - v2.v[2];
  hv[3] = v1.v[3] - v2.v[3];
  return hv;
}

vector4 operator*(double q, const vector4 & v)
{
  vector4 r;
  r[0] = q*v.v[0];
  r[1] = q*v.v[1];
  r[2] = q*v.v[2];
  r[3] = q*v.v[3];
  return r;
}

vector4 operator*(const vector4 & v , double q)
{
  vector4 r;
  r[0] = q*v.v[0];
  r[1] = q*v.v[1];
  r[2] = q*v.v[2];
  r[3] = q*v.v[3];
  return r;
}

vector4 operator/(const vector4 & v , double q)
{
  vector4 r;
  r[0] = v.v[0]/q;
  r[1] = v.v[1]/q;
  r[2] = v.v[2]/q;
  r[3] = v.v[3]/q;
  return r;
}

vector4 operator*(const vector4 & v, const tensor4 & t)
{
  vector4 r;
  r.v[0] = v.v[0]*t.t[0][0] - v.v[1]*t.t[1][0] - v.v[2]*t.t[2][0] - v.v[3]*t.t[3][0];
  r.v[1] = v.v[0]*t.t[0][1] - v.v[1]*t.t[1][1] - v.v[2]*t.t[2][1] - v.v[3]*t.t[3][1];
  r.v[2] = v.v[0]*t.t[0][2] - v.v[1]*t.t[1][2] - v.v[2]*t.t[2][2] - v.v[3]*t.t[3][2];
  r.v[3] = v.v[0]*t.t[0][3] - v.v[1]*t.t[1][3] - v.v[2]*t.t[2][3] - v.v[3]*t.t[3][3];
  return r;
}

vector4 operator*(const tensor4 & t, const vector4 & v)
{
  vector4 r;
  r.v[0] = t.t[0][0]*v.v[0] - t.t[0][1]*v.v[1] - t.t[0][2]*v.v[2] - t.t[0][3]*v.v[3];
  r.v[1] = t.t[1][0]*v.v[0] - t.t[1][1]*v.v[1] - t.t[1][2]*v.v[2] - t.t[1][3]*v.v[3];
  r.v[2] = t.t[2][0]*v.v[0] - t.t[2][1]*v.v[1] - t.t[2][2]*v.v[2] - t.t[2][3]*v.v[3];
  r.v[3] = t.t[3][0]*v.v[0] - t.t[3][1]*v.v[1] - t.t[3][2]*v.v[2] - t.t[3][3]*v.v[3];
  return r;
}

double operator*(const vector4 & v1, const vector4 & v2)
{
  return (v1.v[0]*v2.v[0]-v1.v[1]*v2.v[1]-v1.v[2]*v2.v[2]-v1.v[3]*v2.v[3]);
}

tensor4::tensor4()
{
  int i,j;
  for(i=0;i<4;i++) t[i] = new double [4];
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j] = 0.;
    }
  }
}

tensor4::tensor4(double a00, double a01, double a02, double a03, double a10, double a11, double a12, double a13, double a20, double a21, double a22, double a23, double a30, double a31, double a32, double a33)
{
  int i;
  for(i=0;i<4;i++) t[i] = new double [4];
  t[0][0] = a00;
  t[0][1] = a01;
  t[0][2] = a02;
  t[0][3] = a03;
  t[1][0] = a10;
  t[1][1] = a11;
  t[1][2] = a12;
  t[1][3] = a13;
  t[2][0] = a20;
  t[2][1] = a21;
  t[2][2] = a22;
  t[2][3] = a23;
  t[3][0] = a30;
  t[3][1] = a31;
  t[3][2] = a32;
  t[3][3] = a33;
}

tensor4::tensor4(double a[4][4])
{
  int i,j;
  for(i=0;i<4;i++) t[i] = new double [4];
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j] = a[i][j];
    }
  }
}

tensor4::tensor4(const tensor4 & rhs)
{
  int i,j;
  for(i=0;i<4;i++) t[i] = new double [4];
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j] = rhs.t[i][j];
    }
  }
}

tensor4::~tensor4()
{ 
  int i;
  for(i=0;i<4;i++) delete t[i];
}

tensor4 tensor4::operator=(const tensor4 &rhs)
{
  int i,j;
  if (this != &rhs)
    {
      for (i=0;i<4;i++) delete t[i];
      for (i=0;i<4;i++) t[i] = new double [4];
      for (i=0;i<4;i++)
	for (j=0;j<4;j++)
	  t[i][j] = rhs.t[i][j];
    }
  return *this;
}

double *&tensor4::operator[](int i1)
{   
  if ((i1<0)||(i1>=4)) {
    cout << "\n\a ****  First index of the tensor out of range.  **** \n\n";
    return t[0];
  }
  return t[i1];
}










tensor4 operator*(const tensor4 & a, const tensor4 & b)
{
  tensor4 c;
  int i,j;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      c.t[i][j] = a.t[i][0]*b.t[0][j] - a.t[i][1]*b.t[1][j] - a.t[i][2]*b.t[2][j] - a.t[i][3]*b.t[3][j];
  return c;
}

tensor4 EucProd(const tensor4 & a, const tensor4 & b)
{
  tensor4 c;
  int i,j;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      c.t[i][j] = a.t[i][0]*b.t[0][j] + a.t[i][1]*b.t[1][j] + a.t[i][2]*b.t[2][j] + a.t[i][3]*b.t[3][j];
  return c;
}

vector4 EucProd(const vector4 &v, const tensor4 &t)
{
  vector4 r;
  r.v[0] = v.v[0]*t.t[0][0] + v.v[1]*t.t[1][0] + v.v[2]*t.t[2][0] + v.v[3]*t.t[3][0];
  r.v[1] = v.v[0]*t.t[0][1] + v.v[1]*t.t[1][1] + v.v[2]*t.t[2][1] + v.v[3]*t.t[3][1];
  r.v[2] = v.v[0]*t.t[0][2] + v.v[1]*t.t[1][2] + v.v[2]*t.t[2][2] + v.v[3]*t.t[3][2];
  r.v[3] = v.v[0]*t.t[0][3] + v.v[1]*t.t[1][3] + v.v[2]*t.t[2][3] + v.v[3]*t.t[3][3];
  return r;
}

vector4 EucProd(const tensor4 &t, const vector4 &v)
{
  vector4 r;
  r.v[0] = t.t[0][0]*v.v[0] + t.t[0][1]*v.v[1] + t.t[0][2]*v.v[2] + t.t[0][3]*v.v[3];
  r.v[1] = t.t[1][0]*v.v[0] + t.t[1][1]*v.v[1] + t.t[1][2]*v.v[2] + t.t[1][3]*v.v[3];
  r.v[2] = t.t[2][0]*v.v[0] + t.t[2][1]*v.v[1] + t.t[2][2]*v.v[2] + t.t[2][3]*v.v[3];
  r.v[3] = t.t[3][0]*v.v[0] + t.t[3][1]*v.v[1] + t.t[3][2]*v.v[2] + t.t[3][3]*v.v[3];
  return r;
}



tensor4 diag(const vector4  a)
{
  int i,j;
  tensor4 t;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      t.t[i][j] = 0.;
  for (i=0;i<4;i++) t.t[i][i] = a.v[i];
  return t;
}

tensor4 diag(double a[4])
{
  int i,j;
  tensor4 t;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      t.t[i][j] = 0.;
  for (i=0;i<4;i++) t.t[i][i] = a[i];
  return t;
}

tensor4 diag(double t00, double t11, double t22, double t33)
{
  int i,j;
  tensor4 t;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      t.t[i][j] = 0.;
  t.t[0][0] = t00;
  t.t[1][1] = t11;
  t.t[2][2] = t22;
  t.t[3][3] = t33;
  return t;
}

/* *************************************************************** */










vector3::vector3()
{
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;
}

vector3::vector3(double a1, double a2, double a3)
{
  v[0] = a1;
  v[1] = a2;
  v[2] = a3;
}

vector3::vector3(double a[3])
{
  v[0] = a[0];
  v[1] = a[1];
  v[2] = a[2];
}

vector3::vector3(const vector3 & rhs)
{
  v[0] = rhs.v[0];
  v[1] = rhs.v[1];
  v[2] = rhs.v[2];
}

vector3::~vector3()
{ }

double  &vector3::operator[](int index)
{
  if ((index<1)||(index>=4)){
    cout << "\n\a ****  3-vector index under/overflow  ****\n\n";
    return v[0];
  }
  return v[index-1];
}

vector3 operator+(const vector3 & v1, const vector3 & v2)
{
  vector3 hv;
  hv.v[0] = v1.v[0] + v2.v[0];
  hv.v[1] = v1.v[1] + v2.v[1];
  hv.v[2] = v1.v[2] + v2.v[2];
  return hv;
}

vector3 operator-(const vector3 & v1, const vector3 & v2)
{
  vector3 hv;
  hv.v[0] = v1.v[0] - v2.v[0];
  hv.v[1] = v1.v[1] - v2.v[1];
  hv.v[2] = v1.v[2] - v2.v[2];
  return hv;
}

vector3 operator*(double q, const vector3 & v)
{
  vector3 r;
  r.v[0] = q*v.v[0];
  r.v[1] = q*v.v[1];
  r.v[2] = q*v.v[2];
  return r;
}

vector3 operator*(const vector3 & v , double q)
{
  vector3 r;
  r.v[0] = q*v.v[0];
  r.v[1] = q*v.v[1];
  r.v[2] = q*v.v[2];
  return r;
}


vector3 operator/(const vector3 & v , double q)
{
  vector3 r;
  r.v[0] = v.v[0]/q;
  r.v[1] = v.v[1]/q;
  r.v[2] = v.v[2]/q;
  return r;
}




vector3 operator*(const vector3 & v, const tensor3 & t)
{
  vector3 r;
  r.v[0] = v.v[0]*t.t[0][0] + v.v[1]*t.t[1][0] + v.v[2]*t.t[2][0];
  r.v[1] = v.v[0]*t.t[0][1] + v.v[1]*t.t[1][1] + v.v[2]*t.t[2][1];
  r.v[2] = v.v[0]*t.t[0][2] + v.v[1]*t.t[1][2] + v.v[2]*t.t[2][2];
  return r;
}

vector3 operator*(const tensor3 & t, const vector3 & v)
{
  vector3 r;
  r.v[0] = t.t[0][0]*v.v[0] + t.t[0][1]*v.v[1] + t.t[0][2]*v.v[2];
  r.v[1] = t.t[1][0]*v.v[0] + t.t[1][1]*v.v[1] + t.t[1][2]*v.v[2];
  r.v[2] = t.t[2][0]*v.v[0] + t.t[2][1]*v.v[1] + t.t[2][2]*v.v[2];
  return r;
}

double operator*(const vector3 & v1, const vector3 & v2)
{
  return (v1.v[0]*v2.v[0]+v1.v[1]*v2.v[1]+v1.v[2]*v2.v[2]);
}

tensor3::tensor3()
{
  int i,j;
  for(i=0;i<3;i++) t[i] = new double [3];
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      t[i][j] = 0.;
    }
  }
}

tensor3::tensor3(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22)
{
  int i;
  for(i=0;i<3;i++) t[i] = new double [3];
  t[0][0] = a00;
  t[0][1] = a01;
  t[0][2] = a02;
  t[1][0] = a10;
  t[1][1] = a11;
  t[1][2] = a12;
  t[2][0] = a20;
  t[2][1] = a21;
  t[2][2] = a22;
}

tensor3::tensor3(double a[3][3])
{
  int i,j;
  for(i=0;i<3;i++) t[i] = new double [3];
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      t[i][j] = a[i][j];
    }
  }
}

tensor3::tensor3(const tensor3 & rhs)
{
  int i,j;
  for(i=0;i<3;i++) t[i] = new double [3];
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      t[i][j] = rhs.t[i][j];
    }
  }
}

tensor3::~tensor3()
{ 
  int i;
  for(i=0;i<3;i++) delete t[i];
}

tensor3 tensor3::operator=(const tensor3 &rhs)
{
  int i,j;
  if (this != &rhs)
    {
      for (i=0;i<3;i++) delete t[i];
      for (i=0;i<3;i++) t[i] = new double [3];
      for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	  t[i][j] = rhs.t[i][j];
    }
  return *this;
}

double *&tensor3::operator[](int i1)
{   
  if ((i1<1)||(i1>3)) {
    cout << "\n\a ****  First index of the tensor out of range.  **** \n\n";
    return t[0];
  }
  return t[i1-1];
}


vector3 rotate(tensor3 R,vector3 v)
{
  return (R*v);
}

vector3 rotate(vector3 v, tensor3 R)
{
  return (R*v);
}


tensor3 EulerRotation(double alpha, double beta, double gamma)
{
  tensor3 L;
  
  L.t[0][0] = cos(alpha)*cos(gamma) - sin(alpha)*cos(beta)*sin(gamma);
  L.t[0][1] = sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*sin(gamma);
  L.t[0][2] = sin(beta)*sin(gamma);

  L.t[1][0] = -cos(alpha)*sin(gamma) - sin(alpha)*cos(beta)*cos(gamma);
  L.t[1][1] = - sin(alpha)*sin(gamma) + cos(alpha)*cos(beta)*cos(gamma);
  L.t[1][2] = sin(beta)*cos(gamma);
  
  L.t[2][0] = sin(alpha)*sin(beta);
  L.t[2][1] = -cos(alpha)*sin(beta);
  L.t[2][2] = cos(beta);

  return (L);
}





/* ****************************************************************** */


tensor4 BoostMatrix(vector4 uu )
{
  tensor4 L;
  double v, ut, theta, phi;
  int i;
  vector4 u;

  u[0] = uu[0];
  u[1] = -uu[1];
  u[2] = -uu[2];
  u[3] = -uu[3];

  /* test the vector */
  if ( (fabs((u*u)-1.))>1.e-6 ) {
    cout << " ****  Input vector to boost matrix is not four-velocity. **** \n";
    cout << " ****  u*u = " << u*u << " ***\n";
    return L;
  }

  v = sqrt(u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);

  if (v==0.) 
  {

	L.t[0][0] = L.t[1][1] = L.t[2][2] = L.t[3][3] = 1.;
	L.t[0][1] = L.t[0][2] = L.t[0][3] = 0.;
	L.t[1][0] = L.t[1][2] = L.t[1][3] = 0.;
	L.t[2][0] = L.t[2][1] = L.t[2][3] = 0.;
	L.t[3][0] = L.t[3][1] = L.t[3][2] = 0.;
  }
  else
  {
  theta = acos(u[3]/v);
	  
  if ( theta<1.e-5 ) 
    {
      phi = 0.;
    }
  else
    {
      ut = sqrt(u[1]*u[1]+u[2]*u[2]);
		if (ut == 0)
			phi = 0.;
		else 
		{
			if (u[2] > 0)
			{
				phi = acos(u[1]/ut);
			}
			else 
			{
				phi = 2.*3.1415926-acos(u[1]/ut);
			}
		}
    }
		
  
  L.t[0][0] = u[0];
  for (i=1;i<4;i++) L.t[i][0] = L.t[0][i] = -u[i];
 
  L.t[1][1] = cos(theta)*cos(theta)*cos(phi)*cos(phi) 
    + u[0]*cos(phi)*cos(phi)*sin(theta)*sin(theta) + sin(phi)*sin(phi);
 
  L.t[1][2] = L.t[2][1] = (u[0]-1.)*cos(phi)*sin(phi)*sin(theta)*sin(theta);

  L.t[1][3] = L.t[3][1] = (u[0]-1.)*cos(theta)*sin(theta)*cos(phi);
  
  L.t[2][2] = cos(phi)*cos(phi) 
    + (cos(theta)*cos(theta) + u[0]*sin(theta)*sin(theta))*sin(phi)*sin(phi);

  L.t[2][3] = L.t[3][2] = (u[0]-1.)*cos(theta)*sin(theta)*sin(phi);

  L.t[3][3] = u[0]*cos(theta)*cos(theta) + sin(theta)*sin(theta);
  }


  return L;
    
}


vector4 boost(vector4 v, tensor4 L)
{
  vector4 boostedv;
  double invariant, dev; 

  invariant = v*v;
  boostedv = EucProd(v,L);

  if ((dev = fabs((boostedv*boostedv - invariant)/invariant)) > 1.e-8)
    {
      //      cerr << " *** WARNING BOOST: fixing 0th component; deviation was "<< dev<<" ***\n";
      boostedv[0] = sqrt(invariant + boostedv[1]*boostedv[1] 
			 + boostedv[2]*boostedv[2] + boostedv[3]*boostedv[3]);
    }
  return boostedv;
}

tensor4 boost(tensor4 a, tensor4 L)
{
  tensor4 b;
  b = EucProd(L,EucProd(a,L));
  return b;
}


void prn(tensor4 L)
{
  cout << "(";
  cout.width(12);
  cout << L[0][0]; 
  cout.width(12);
  cout << L[0][1]; 
  cout.width(12);
  cout << L[0][2]; 
  cout.width(12);
  cout << L[0][3]; 
  cout << " ) \n";
  cout << "( ";
  cout.width(12);
  cout << L[1][0]; 
  cout.width(12);
  cout << L[1][1]; 
  cout.width(12);
  cout << L[1][2]; 
  cout.width(12);
  cout << L[1][3]; 
  cout << " ) \n";
  cout << "( ";
  cout.width(12);
  cout << L[2][0]; 
  cout.width(12);
  cout << L[2][1]; 
  cout.width(12);
  cout << L[2][2]; 
  cout.width(12);
  cout << L[2][3]; 
  cout << " ) \n";
  cout << "( ";
  cout.width(12);
  cout << L[3][0]; 
  cout.width(12);
  cout << L[3][1]; 
  cout.width(12);
  cout << L[3][2]; 
  cout.width(12);
  cout << L[3][3]; 
  cout << " ) \n";
}


void prn(vector4 a)
{
  cout << "("; 
  cout.width(12);
  cout << a[0];
  cout << ",  ",
  cout.width(12);
  cout << a[1];
  cout << ",  ",
  cout.width(12);
  cout << a[2];
  cout << ",  ",
  cout.width(12);
  cout << a[3];
  cout << ")\n";
}

void prt(vector4 a)
{
  cout << "("; 
  cout.width(12);
  cout << a[0];
  cout << ",  ",
  cout.width(12);
  cout << a[1];
  cout << ",  ",
  cout.width(12);
  cout << a[2];
  cout << ",  ",
  cout.width(12);
  cout << a[3];
  cout << ")";
}


/* *** define the class particle *** */

ParticleReg::ParticleReg()
{
  x[0] = x[1] = x[2] = x[3] = 0.;
  p[0] = p[1] = p[2] = p[3] = 0.;
  id = 0;
  droplet = -1;
  resdecay = 0;
}

ParticleReg::ParticleReg(const ParticleReg & rhs)
{
  x = rhs.x;
  p = rhs.p;
  id = rhs.id;
  droplet = rhs.droplet;
  decay = rhs.decay;
  resdecay = rhs.resdecay;
}

ParticleReg::ParticleReg(int iden, vector4 lsp, vector4 pi)
{
  x[0] = lsp[0];
  x[1] = lsp[1];
  x[2] = lsp[2];
  x[3] = lsp[3];
  p[0] = pi[0];
  p[1] = pi[1];
  p[2] = pi[2];
  p[3] = pi[3];
  id = iden;
  droplet = -1;
  resdecay = 0;
}

ParticleReg::ParticleReg(int iden, vector4 lsp, vector4 pi, int drp)
{
  x[0] = lsp[0];
  x[1] = lsp[1];
  x[2] = lsp[2];
  x[3] = lsp[3];
  p[0] = pi[0];
  p[1] = pi[1];
  p[2] = pi[2];
  p[3] = pi[3];
  id = iden;
  droplet = drp;
  resdecay = 0;
}


ParticleReg::ParticleReg(int iden, vector4 lsp, vector4 pi, int drp, bool rd)
{
  x[0] = lsp[0];
  x[1] = lsp[1];
  x[2] = lsp[2];
  x[3] = lsp[3];
  p[0] = pi[0];
  p[1] = pi[1];
  p[2] = pi[2];
  p[3] = pi[3];
  id = iden;
  droplet = drp;
  resdecay = rd;
}


ParticleReg::~ParticleReg()
{ }


double ParticleReg::mass() const
{
  return sqrt(p*p);
}


/* ***************************************** */


double SR_redI0(double x)

/*  this computes the Bessel function I0 multiplied with exp(-x)  */

{
  double ax, ans, y;

  if ((ax=fabs(x))<3.75){
    y = x/3.75;
    y *= y;

    ans = 1.+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(.2659732+y*(.0360768+
               y*.45813e-2)))));
    ans *= exp(-ax);
  } else {
    y = 3.75/x;
    ans = (1./sqrt(ax))*
      (.39894228+y*(.01328592+y*(.225319e-2+y*
      (-0.157565e-2+y*(0.916281e-2+y*(-0.02057706+y*
      (0.02635537+y*(-0.01647633+y*0.392377e-2))))))));
  }

  return ans;
}


double SR_redI1(double x)

/*   This computes the Bessel function I1 multiplied with exp(-x)  */

{
  double ax, ans, y;

  if ((ax=fabs(x))<3.75) {
    y = x/3.75;
    y *= y;
    ans = ax*
      (0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.02658733+y*
      (0.301532e-2+y*0.32411e-3))))));
    ans *= exp(-ax);
  } else {
    y = 3.75/ax;
    ans = (1./sqrt(ax))*
      (0.39894228+y*(-0.03988024+y*(-0.362018e-2+y*(0.163801e-2+y*
      (-0.01031555+y*(0.02282967+y*(-0.02895312+y*
      (0.01787654-y*0.420059e-2))))))));
  }

  return ans;
}



double SR_redK0(double x)

/*   This computes the Bessel function K0 multipied with exp(x)  */

{

  double y, ans;

  if (x <= 2.){
    y = x*x/4.;
    ans = -log(x/2.)*SR_redI0(x)*exp(x) - 0.57721566 + y*
      (0.42278420 + y*(0.23069756+y*(0.03488590+y*(0.262698e-2+y*
                                                   (0.1075e-3+y*0.74e-5)))));
    ans *= exp(x);
  } else {
    y = 2./x;
    ans = (1./sqrt(x))*(1.25331414 + y*
      (-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*
                                        (0.587873e-2+y*
                                         (-0.251540e-2+y*0.53208e-3))))));
  }

  return ans;
}


double SR_redK1(double x)

/*   This computes the Bessel function K1 multiplied with exp(x)  */

{
  //  double reedI1();
  double y, ans;

  if ( x <= 2.){
    y = x*x/4.;
    ans = log(x/2.)*SR_redI1(x)*exp(x) + 1./x + y/x *
      (0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.01919402+y*
      (-0.110404e-2-y*0.4686e-4)))));
    ans *= exp(x);
  }  else {
    y = 2./x;
    ans = (1./sqrt(x))*
      (1.25331414+y*(0.23498619+y*(-0.03655620+y*(0.01504268+y*
      (-0.780353e-2+y*(0.325614e-2-y*0.68245e-3))))));
  }

  return ans;

}

double sum (int index, int op,double *massl){
double s=0;
for (int i=index;i<=op;i++){
s = s+massl[i];
}
return s;
}


double SR_redK2(double x)

/*   This gives the Bessel function K2 multiplied with exp(x)  */

{
  double ans;

  ans = (2./x)*SR_redK1(x) + SR_redK0(x);

  return ans;
}
