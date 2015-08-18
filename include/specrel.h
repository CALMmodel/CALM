/* relativistic algebra */
/* we define all vectors and tensors as covariant                    */
/* Minkowski metric is taken into account in the definitions         */
/* of operation * , which is done by overloading the oparator        */
#ifndef _SPECREL_H
#define _SPECREL_H

#define SRPI 3.141592653589793238

#define SRsign(x) ( (x)>=0 ? 1 : -1 )

class tensor4
{
public:
  tensor4();
  tensor4(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
  tensor4(double a[4][4]);
  tensor4(const tensor4 &);
  ~tensor4();
  tensor4 operator=(const tensor4 &);
  double * t[4];
  double * &operator[](int);
  friend tensor4 operator*(const tensor4 &, const tensor4 &);
private:
}; 


class vector4
{
public:
  vector4();
  vector4(double a0, double a1, double a2, double a3);
  vector4(double a[4]);
  vector4(const vector4 &);
  ~vector4();
  double v[4];
  double &operator[](int);
  friend vector4 operator+(const vector4 &, const vector4 &);
  friend vector4 operator-(const vector4 &, const vector4 &);
  friend vector4 operator*(double , const vector4 &);
  friend vector4 operator*(const vector4 & , double);
  friend vector4 operator/(const vector4 & , double);
  friend vector4 operator*(const vector4 &, const tensor4 &);// row vector x tensor
  friend vector4 operator*(const tensor4 &,const vector4 &);// tensor x column vector
private:
};

/* this must give the scalar procuct of four-vectors in Minkowski metric */
double operator*(const vector4&, const vector4&);

/* these functions return products in Euclidean metrics */
tensor4 EucProd(const tensor4 &, const tensor4 &);
vector4 EucProd(const vector4 &, const tensor4 &);
vector4 EucProd(const tensor4 &, const vector4 &);

/* this produces a diagonal four-tensor */
tensor4 diag(const vector4 );
tensor4 diag(double a[4]);
tensor4 diag(double, double, double, double);

/* this calculates the boost matrix (to be used from left on 4-vectors) from given boost velocity */
tensor4 BoostMatrix(vector4 );
// tensor4 BoostMatrix(double v[3]);
// tensor4 BoostMatrix(double , double , double , double);

vector4 boost(vector4 v, tensor4 L);
tensor4 boost(tensor4 a, tensor4 L);

void prn(tensor4);
void prn(vector4);
void prt(vector4);





//class for three-matrices
class tensor3
{
public:
  tensor3();
  tensor3(double, double, double, double, double, double, double, double, double);
  tensor3(double a[3][3]);
  tensor3(const tensor3 &);
  ~tensor3();
  tensor3 operator=(const tensor3 &);
  double * t[3];
  double * &operator[](int);
  friend tensor3 operator*(const tensor3 &, const tensor3 &);
private:
}; 





// class for three-vectors
class vector3
{
public:
  vector3();
  vector3(double a1, double a2, double a3);
  vector3(double a[3]);
  vector3(const vector3 &);
  ~vector3();
  double v[3];
  double &operator[](int);
  friend vector3 operator+(const vector3 &, const vector3 &);
  friend vector3 operator-(const vector3 &, const vector3 &);
  friend vector3 operator*(double , const vector3 &);
  friend vector3 operator*(const vector3 & , double);
  friend vector3 operator/(const vector3 & , double);
  friend vector3 operator*(const vector3 &, const tensor3 &);// row vector x tensor
  friend vector3 operator*(const tensor3 &,const vector3 &);// tensor x column vector
private:
};



vector3 rotate(tensor3 R,vector3 v);
vector3 rotate(vector3 v, tensor3 R);
tensor3 EulerRotation(double alpha, double beta, double gamma);


// class for defining particles
class ParticleReg
{
public:
  ParticleReg();
  ParticleReg(int, vector4 , vector4);
  ParticleReg(int, vector4 , vector4, int);
  ParticleReg(int, vector4 , vector4, int, bool);
  ParticleReg(const ParticleReg & rhs);
  ~ParticleReg();
  vector4 x;
  vector4 p;
  long int id;
  int droplet; // arbitrary, number of the emitting droplet can be stored here
  bool resdecay; // 1 if coming from resonance decay, 0 otherwise 
  int decay;
  double mass() const;
private:
};

// reduced Bessel functions
double SR_redI0(double);
double SR_redI1(double);
double SR_redK0(double);
double SR_redK1(double);
double SR_redK2(double);
double sum (int index, int op, double *massl);

#endif
