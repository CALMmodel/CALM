#include "reggae.h"

// number of collisions per particle in the rescattering phase 
// divided by 2
int RG_opak = 4;

using namespace std;
/* the random generator routine */

double KAS_rndm(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[KAS_RND_NTAB];
  double temp;

  if (*idum <= 0){
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j=KAS_RND_NTAB+7;j>=0;j--){
      k=(*idum)/KAS_RND_IQ1;
      *idum = KAS_RND_IA1*(*idum-k*KAS_RND_IQ1)-k*KAS_RND_IR1;
      if (*idum < 0) *idum += KAS_RND_IM1;
      if (j < KAS_RND_NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/KAS_RND_IQ1;
  *idum = KAS_RND_IA1*(*idum-k*KAS_RND_IQ1)-k*KAS_RND_IR1;
  if (*idum < 0) *idum += KAS_RND_IM1;
  k = idum2/KAS_RND_IQ2;
  idum2 = KAS_RND_IA2*(idum2-k*KAS_RND_IQ2)-k*KAS_RND_IR2;
  if (idum2 < 0) idum2 += KAS_RND_IM2;
  j = iy/KAS_RND_NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += KAS_RND_IMM1;
  if ((temp=KAS_RND_AM*iy) > KAS_RND_RNMX) return KAS_RND_RNMX;
  else return temp;
}


void quicksort(int n, double a[])
{
  int left, right, middle, length;
  int part;
  int i,j;
  int pending[100];
  int onstack;
  double temp;

  pending[0] = 0;
  pending[1] = n-1;
  onstack = 1;

  while (onstack>=0)
    {
      right = pending[onstack--];
      left = pending[onstack--];
      length = right - left + 1;

      if (length < SortLength)
	{// do the straight insertion
	  for (i=left+1;i<=right;i++)
	    {
	      temp = a[i];
	      j = i-1;
	      while ( (temp < a[j]) && (j >= left))
		{
		  a[j+1] = a[j];
		  j -= 1;
		}
	      a[j+1] = temp;
	    }
	}
      else
	{// if partition is too large, do real quicksort
	  middle = left + 1;
	  // now pick the median of three values
	  if (a[middle]<a[left]) 
	    {temp = a[middle]; a[middle]=a[left]; a[left]=temp;}
	  if (a[right]<a[middle])
	    {temp = a[right]; a[right] = a[middle]; a[middle] = temp;}
	  if (a[middle]<a[left]) 
	    {temp = a[middle]; a[middle]=a[left]; a[left]=temp;}
	  // so we pick a[middle]

	  i = left+1;
	  j = right;

	  for (;;)
	    {
	      do i++; while (a[i]<a[middle]);
	      do j--; while (a[j]>a[middle]);
	      if (i>j) break;
	      // switch the entries
	      temp = a[i]; a[i] = a[j]; a[j] = temp;
	    }
	  //insert partitioning element
	  temp = a[middle]; a[middle] = a[j]; a[j] = temp;
	  part = j;
	  // stack new partitions
	  onstack += 1;
	  pending[onstack] = part+1;
	  onstack += 1;
	  pending[onstack] = right;
	  onstack += 1;
	  pending[onstack] = left;
	  onstack += 1;
	  pending[onstack] = part-1;
	}
    }

}

double Mconserv (vector4 P, int n,double *mass, vector4 *op, long int *seed){
  double *arrnd = (double*)malloc(n * sizeof (double));
  double *aimass = (double*)malloc(n * sizeof (double));
  int i,l,j;
  double M,ps,w,sum1;
  double rnd,theta;
  double oarg;

  vector4 pa,Q,pdec,sku;
  vector4 bv;
  tensor4 bm;


  for (i=0;i<(n-2);i++){
    arrnd[i]=KAS_rndm(seed);
  }
  quicksort(n-2,arrnd);
  M=sqrt(P*P);	//invariant mass
  sum1 = M - sum(0,n,mass);
  for (int i=1;i<n-1;i++){
    aimass[i] = (arrnd[i-1]*(sum1)) + sum(0,i,mass);
  }
  aimass[n-1]=M;
  aimass[0]=mass[0];


  pdec=P;
  j=n-1;
  while (j>0) {

    oarg = (aimass[j]*aimass[j] - (aimass[j-1]*aimass[j-1] + mass[j]*mass[j]))
      * (aimass[j]*aimass[j] - (aimass[j-1]*aimass[j-1] + mass[j]*mass[j]))
      - 4.*mass[j]*mass[j]*aimass[j-1]*aimass[j-1];
	
    if (-1.e-6<oarg && oarg<0.) oarg = 0.;
	
    ps= sqrt(oarg)
      / (2.*aimass[j]); 
	
    // generate angles
    rnd = 2.*DGPI*KAS_rndm(seed);
    theta = acos(1.-2.*KAS_rndm(seed));

    // calculate momenta
    pa[1] = ps * sin(theta) * cos(rnd);
    pa[2] = ps * sin(theta) * sin(rnd);
    pa[3] = ps * cos(theta);
    pa[0] = sqrt(mass[j]*mass[j] 
		 + ps*ps);

    Q[1] = -pa[1];
    Q[2] = -pa[2];
    Q[3] = -pa[3];
    Q[0] = sqrt(aimass[j-1]*aimass[j-1] 
		+ ps*ps);

    //now moment must be boosted
    //get boost-velocity
    bv=pdec/aimass[j];

    //calculate boost-tensor
    bm=BoostMatrix(bv);


    //boost the momenta
    pa=boost(pa,bm);
    op[j]=pa;
    Q=boost(Q,bm);
    pdec=Q;

    j--;

 
  }
  op[0]=Q;

  w=1.;
  w = (pow(SRPI,n-1))/(2*M);
  for (int l=1;l<n;l++){
    w *= sqrt(aimass[l]*aimass[l] 
	      + pow((aimass[l-1]*aimass[l-1] - mass[l]*mass[l])/aimass[l],2)
	      -2*(aimass[l-1]*aimass[l-1] + mass[l]*mass[l]));
		
  }
  free(arrnd);
  free(aimass);
  return w;
}

double collision (int n,vector4 *avec, long int *seed)
{
  int j;
  double rnd,theta;
  double psp1,psp2,sqs,r,r2;
  vector4 p1,p2,ptot,vboost,p1star,p2star;
  tensor4 lambda;
  double m12, m22;
  double scms;
  for (int c=0;c<RG_opak;c++){
    int ran = static_cast<int>(n*KAS_rndm(seed));
    while (ran == 0 ) {
      ran = static_cast<int>(n*KAS_rndm(seed));}
    //cout << "posun je " << ran << "\n"; 
    for (int i=0;i<n;i++){
      j=i+ran;
      if (j>n-1) j=j-n;
      // cout << i <<" " << j << "\n";
			
      //tu zacina zrazka
      p1=avec[i];
      p2=avec[j];
      m12 = p1*p1;
      m22 = p2*p2;
      ptot=p1+p2;
      scms = ptot*ptot;
      // collide if s is not too small 
      if ((scms - (m12+m22+2.*sqrt(m12*m22)))/(p1[0]+p2[0])/(p1[0]+p2[0]) > 0.0001) 
	{
	  sqs=sqrt(scms);
	  // cout << scms << " I am here \n";
	  vboost=ptot/sqs;

				vboost[1]=-vboost[1];
				vboost[2]=-vboost[2];
				vboost[3]=-vboost[3];
				lambda=BoostMatrix(vboost);
				p1star=boost(p1,lambda);
				p2star=boost(p2,lambda);
 
				// generate angles
				rnd = 2.*DGPI*KAS_rndm(seed);
				theta = acos(1.-2.*KAS_rndm(seed));
		
				r=sqrt(p1star[1]*p1star[1]+p1star[2]*p1star[2]+p1star[3]*p1star[3]);
				r2=sqrt(p2star[1]*p2star[1]+p2star[2]*p2star[2]+p2star[3]*p2star[3]);

				// calculate momenta
      			p1star[1] = r * sin(theta) * cos(rnd);
      			p1star[2] = r * sin(theta) * sin(rnd);
      			p1star[3] = r * cos(theta);
      		
				p2star[1] = -p1star[1];
				p2star[2] = -p1star[2];
				p2star[3] = -p1star[3];

				vboost[1]=-vboost[1];
				vboost[2]=-vboost[2];
				vboost[3]=-vboost[3];
				lambda=BoostMatrix(vboost);
			
				avec[i]=boost(p1star,lambda);
				avec[j]=boost(p2star,lambda);
				// restore masses
				avec[i][0] = sqrt(avec[i][1]*avec[i][1]+avec[i][2]*avec[i][2]
								  +avec[i][3]*avec[i][3]+m12);
				avec[j][0] = sqrt(avec[j][1]*avec[j][1]+avec[j][2]*avec[j][2]
								  +avec[j][3]*avec[j][3]+m22);
			}
		}
	}
}
