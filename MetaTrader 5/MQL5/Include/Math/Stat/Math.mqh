//+------------------------------------------------------------------+
//|                                                         Math.mqh |
//|                        Copyright 2016, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2016, MetaQuotes Software Corp."
#property link      "https://www.mql5.com"

#define ERR_OK                0
#define ERR_ARGUMENTS_NAN     1
#define ERR_ARGUMENTS_INVALID 2
#define ERR_RESULT_INFINITE   3
#define ERR_NON_CONVERGENCE   4
//+------------------------------------------------------------------+
//| Nan                                                              |
//+------------------------------------------------------------------+
double Nan(long bit_value)
  {
   struct L { long   x; } l; l.x=bit_value;
   struct D { double x; } d=(D)l;
   return(d.x);
  }

double QNaN   =Nan(0x7FF7000000000000);   // QNaN 
double QPOSINF=Nan(0x7FF0000000000000);   // +infinity 1.#INF 
double QNEGINF=Nan(0xFFF0000000000000);   // -infinity -1.#INF 
//+------------------------------------------------------------------+
//| MathRandom                                                       |
//+------------------------------------------------------------------+
double MathRandomNonZero(void)
  {
   double rnd=0;
//--- except 0 and 1
   while(rnd==0.0 || rnd==1.0)
      rnd=MathRand()/32767.0;
//---
   return(rnd);
  }
//+------------------------------------------------------------------+
//| Computes the minimum value in array[]                            |
//+------------------------------------------------------------------+
double MathMin(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- first element by default, find minimum
   double min_value=array[ind1];

   for(int i=ind1+1; i<=ind2; i++)
      min_value=MathMin(min_value,array[i]);
//--- return minimum value
   return(min_value);
  }
//+------------------------------------------------------------------+
//| Computes the maximum value in array[]                            |
//+------------------------------------------------------------------+
double MathMax(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- first element by default, find maximum
   double max_value=array[ind1];

   for(int i=ind1+1; i<=ind2; i++)
      max_value=MathMax(max_value,array[i]);
//--- return maximum value
   return(max_value);
  }
//+------------------------------------------------------------------+
//| Computes the range of the values in array[]                      |
//+------------------------------------------------------------------+
double MathRange(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- default values, find minimum and maximum values
   double min_value=array[ind1];
   double max_value=array[ind1];

   for(int i=ind1+1; i<=ind2; i++)
     {
      double value=array[i];
      min_value=MathMin(min_value,value);
      max_value=MathMax(max_value,value);
     }
//--- return range
   return(max_value-min_value);
  }
//+------------------------------------------------------------------+
//| Computes the sum of the values in array[]                        |
//+------------------------------------------------------------------+
double MathSum(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate sum
   double sum=0.0;
   for(int i=ind1; i<=ind2; i++)
      sum+=array[i];
//--- return sum
   return(sum);
  }
//+------------------------------------------------------------------+
//| Computes the standard deviation of the values in array[]         |
//+------------------------------------------------------------------+
double MathStandardDeviation(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<=1)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
//--- average mean
   mean=mean/data_count;
//--- calculate standard deviation   
   double sdev=0;
   for(int i=ind1; i<=ind2; i++)
      sdev+=MathPow(array[i]-mean,2);
//--- return standard deviation
   return MathSqrt(sdev/(data_count-1));
  }
//+------------------------------------------------------------------+
//| Computes the average absolute deviation of the values in array[] |
//+------------------------------------------------------------------+
double MathAverageDeviation(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<=1)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
//--- calculate average deviation
   double adev=0;
   for(int i=ind1; i<=ind2; i++)
      adev+=MathAbs(array[i]-mean);
   adev=adev/data_count;
//--- return average deviation
   return(adev);
  }
//+------------------------------------------------------------------+
//| Computes the median value of the values in array[]               |
//+------------------------------------------------------------------+
double MathMedian(double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(QNaN);
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- prepare sorted values
   double sorted_values[];
   ArrayCopy(sorted_values,array,0,start,count);
   ArraySort(sorted_values);
//--- calculate median for odd and even cases
//--- data_count=odd
   if(data_count%2==1)
      return(sorted_values[data_count/2]);
   else
//--- data_count=even
      return(0.5*(sorted_values[(data_count-1)/2]+sorted_values[(data_count+1)/2]));
  }
//+------------------------------------------------------------------+
//| Computes the mean value of the values in array[]                 |
//+------------------------------------------------------------------+
double MathMean(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<1)
      return(QNaN); // need at least 1 observation
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
//--- return mean
   return(mean);
  }
//+------------------------------------------------------------------+
//| Computes the variance of the values in array[]                   |
//+------------------------------------------------------------------+
double MathVariance(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<2)
      return(QNaN); // need at least 2 observations
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
//--- calculate variance
   double variance=0;
   for(int i=ind1; i<=ind2; i++)
      variance+=MathPow(array[i]-mean,2);
   variance=variance/(data_count-1);
//--- return variance
   return(variance);
  }
//+------------------------------------------------------------------+
//| Computes the skewness of the values in array[]                   |
//+------------------------------------------------------------------+
double MathSkewness(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<3)
      return(QNaN); // need at least 3 observations
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
//--- calculate variance and skewness
   double variance=0;
   double skewness=0;
   for(int i=ind1; i<=ind2; i++)
     {
      double sqr_dev=MathPow(array[i]-mean,2);
      skewness+=sqr_dev*(array[i]-mean);
      variance+=sqr_dev;
     }
   variance=(variance)/(data_count-1);
   double v3=MathPow(MathSqrt(variance),3);
//---
   if(v3!=0)
     {
      skewness=skewness/(data_count*v3);
      //--- return skewness
      return(skewness);
     }
   else
      return(QNaN);
  }
//+------------------------------------------------------------------+
//| Computes the kurtosis of the values in array[]                   |
//+------------------------------------------------------------------+
double MathKurtosis(const double &array[],const int start=0,const int count=WHOLE_ARRAY)
  {
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count<4)
      return(QNaN); // need at least 4 observations
   if(start+data_count>size)
      return(QNaN);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   double mean=0.0;
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
//--- calculate variance and kurtosis
   double variance=0;
   double kurtosis=0;
   for(int i=ind1; i<=ind2; i++)
     {
      double sqr_dev=MathPow(array[i]-mean,2);
      variance+=sqr_dev;
      kurtosis+=sqr_dev*sqr_dev;
     }
//--- calculate variance
   variance=(variance)/(data_count-1);
   double v4=MathPow(MathSqrt(variance),4);

   if(v4!=0)
     {
      //--- calculate kurtosis
      kurtosis=kurtosis/(data_count*v4);
      kurtosis-=3;
      //--- return kurtosis
      return(kurtosis);
     }
   else
      return(QNaN);
  }
//+------------------------------------------------------------------+
//| Computes 4 first moments of the values in array[]                |
//+------------------------------------------------------------------+
bool MathMoments(const double &array[],double &mean,double &variance,double &skewness,double &kurtosis,const int start=0,const int count=WHOLE_ARRAY)
  {
//--- initial values
   mean=0.0;
   variance=0.0;
   skewness=0.0;
   kurtosis=0.0;
//--- get data size
   int size=ArraySize(array);
   int data_count=0;
//--- set data count
   if(count==WHOLE_ARRAY)
      data_count=size;
   else
      data_count=count;
//--- check data range
   if(data_count==0)
      return(false);
   if(start+data_count>size)
      return(false);
//--- set indexes
   int ind1=start;
   int ind2=ind1+data_count-1;
//--- calculate mean
   for(int i=ind1; i<=ind2; i++)
      mean+=array[i];
   mean=mean/data_count;
   double S=0.0;
//--- 
   if(data_count>1)
     {
      //--- calculate variance
      for(int i=ind1; i<=ind2; i++)
         variance+=MathPow(array[i]-mean,2);
      variance=variance/(data_count-1);
      S=MathSqrt(variance);
     }
   else variance=QNaN;
//---
   if(data_count>2 && S>0.0)
     {
      //--- calculate skewness
      for(int i=ind1; i<=ind2; i++)
         skewness+=MathPow(array[i]-mean,3);
      skewness=skewness/(data_count*MathPow(S,3));
     }
   else skewness=QNaN;
//---
   if(data_count>3 && S>0.0)
     {
      //--- calculate kurtosis
      for(int i=ind1; i<=ind2; i++)
         kurtosis+=MathPow(array[i]-mean,4);
      kurtosis=kurtosis/(data_count*MathPow(S,4));
      kurtosis-=3;
     }
   else kurtosis=QNaN;
//--- check values and return calculation result
   if((!MathIsValidNumber(variance)) || (!MathIsValidNumber(skewness)) || (!MathIsValidNumber(kurtosis)))
      return(false);
   else
      return(true);
  }

#define M_1_SQRT_2PI  1/MathSqrt(2*M_PI)

double FactorialsTable[21]=
  {
   1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,
   6227020800,87178291200,1307674368000,20922789888000,355687428096000,
   6402373705728000,121645100408832000,2432902008176640000
  };
//+------------------------------------------------------------------+
//| MathPowInt                                                       |
//+------------------------------------------------------------------+
double MathPowInt(const double x,const int power)
  {
//--- check power
   if(power==0)
      return(1.0);
   if(power<0)
      return(0);
//--- calculate product
   double val=x;
   for(int i=2; i<=power; i++)
      val*=x;
//--- return result
   return(val);
  }
//+------------------------------------------------------------------+
//| MathFactorial                                                    |
//+------------------------------------------------------------------+
double MathFactorial(const int n)
  {
   if(n<0)
      return(0);
   if(n<=20)
      //--- use values from the factorials table
      return(FactorialsTable[n]);
   else
     {
      //--- calculate product starting from 20th element of factorials table
      double val=FactorialsTable[20];
      for(int i=21; i<=n; i++)
         val*=i;
      //--- return result
      return(val);
     }
  }
//+------------------------------------------------------------------+
//| MathTrunc                                                        |
//+------------------------------------------------------------------+
double MathTrunc(const double x)
  {
   if(x>=0)
      return(MathFloor(x));
   else
      return(MathCeil(x));
  }
//+------------------------------------------------------------------+
//| MathRound                                                        |
//+------------------------------------------------------------------+
//| Round a double number to a given precision                       |
//+------------------------------------------------------------------+
double MathRound(const double x,const int digits)
  {
   if(!MathIsValidNumber(x))
      return(QNaN);

   if(x==0.0)
      return(x);

   if(digits>0)
     {
      double sign=1.0;
      double xx=x;
      if(xx<0.0)
        {
         xx=-xx;
         sign=-1.0;
        }
      double pwr_factor=MathPow(10,digits);
      double int_value=MathFloor(xx);
      double res=(xx-int_value)*pwr_factor;
      res=MathFloor(MathAbs(res)+0.5);
      return(sign*(int_value+res/pwr_factor));
     }
   else
   if(digits==0)
      return MathRound(x);
//---
   return(0);
  }
//+------------------------------------------------------------------+
//| MathArctan2                                                      |
//+------------------------------------------------------------------+
//| The function calculates polar angle of vector (x,y)              |
//| in range [-pi, pi]. It uses signs to determine quadrant.         |
//+------------------------------------------------------------------+
double MathArctan2(const double y,const double x)
  {
   if(x>0.0)
      return MathArctan(y/x);
   else
   if(x<0.0)
     {
      if(y>0.0)
         return(MathArctan(y/x)+M_PI);
      else
         return(MathArctan(y/x)-M_PI);
     }
//--- case x==0.0
   if(y>0.0)
      return(M_PI_2);
   else
   if(y<0.0)
      return(-M_PI_2);
   else
      return(QNaN);
  }
//+------------------------------------------------------------------+
//| Comments from original FORTRAN code:                             |
//| http://www.netlib.org/specfun/gamma                              |
//|                                                                  |
//| This routine calculates the GAMMA function for a real argument X.|
//| Computation is based on an algorithm outlined in reference 1.    |
//| The program uses rational functions that approximate the GAMMA   |
//| function to at least 20 significant decimal digits. Coefficients |
//| for the approximation over the interval (1,2) are unpublished.   |
//| Those for the approximation for X .GE. 12 are from reference 2.  |
//|                                                                  |
//| References:                                                      |
//| 1. "An Overview of Software Development for Special Functions"   |
//|    W. J. Cody, Lecture Notes in Mathematics, 506,                |
//|    Numerical Analysis Dundee, 1975, G. A. Watson (ed.)           |
//|    Springer Verlag, Berlin, 1976.                                |
//| 2. Computer Approximations, Hart, Et. Al., Wiley and sons,       |
//|    New York, 1968.                                               |
//|                                                                  |
//| Authors:                                                         |
//|    W. J. Cody and L. Stoltz, Applied Mathematics Division        |
//|    Argonne National Laboratory, Argonne, IL 60439                |
//+------------------------------------------------------------------+
double MathGamma(const double x)
  {
//--- mathematical constants
   const double logsqrt2pi=MathLog(MathSqrt(2*M_PI));
//--- machine dependent parameters
   const double xminin=2.23e-308;
   const double xbig=171.624;
   double eps=2.22e-16;
//--- numerator and denominator coefficients for rational minimax
//--- approximation over (1,2)
   const double p[8]=
     {
      -1.71618513886549492533811e+0,2.47656508055759199108314e+1,
      -3.79804256470945635097577e+2,6.29331155312818442661052e+2,
      8.66966202790413211295064e+2,-3.14512729688483675254357e+4,
      -3.61444134186911729807069e+4,6.64561438202405440627855e+4
     };
   const double q[8]=
     {
      -3.08402300119738975254353e+1,3.15350626979604161529144e+2,
      -1.01515636749021914166146e+3,-3.10777167157231109440444e+3,
      2.25381184209801510330112e+4,4.75584627752788110767815e+3,
      -1.34659959864969306392456e+5,-1.15132259675553483497211e+5
     };
//--- coefficients for minimax approximation over (12, inf)
   const double c[7]=
     {
      -1.910444077728e-03,8.4171387781295e-04,
      -5.952379913043012e-04,7.93650793500350248e-04,
      -2.777777777777681622553e-03,8.333333333333333331554247e-02,
      5.7083835261e-03
     };

   int    parity=0;
   double fact=1.0;
   int    n=0;
   double y=x,y1,res,z,xnum=0,xden=0,ysq=0,sum=0;

   if(y<=0.0)
     {
      //--- argument is negative
      y=-x;
      y1=MathTrunc(y);
      res=y-y1;
      if(res!=0.0)
        {
         if(y1!=MathTrunc(y1*0.5)*2.0)
            parity=1;
         fact=-M_PI/MathSin(M_PI*res);
         y=y+1.0;
        }
      else
        {
         return(QPOSINF);
        }
     }
//--- argument is positive
   if(y<eps)
     {
      //--- argument < eps
      if(y>=xminin)
        {
         res=1.0/y;
        }
      else
        {
         return(QPOSINF);
        }
     }
   else
   if(y<12.0)
     {
      y1=y;
      if(y<1.0)
        {
         //--- 0.0 < argument < 1.0
         z = y;
         y = y + 1.0;
        }
      else
        {
         //--- 1.0 < argument < 12.0, reduce argument if necessary
         n = int(y) - 1;
         y = y - double(n);
         z = y - 1.0;
        }
      //--- evaluate approximation for 1.0 < argument < 2.0
      xnum=0.0;
      xden=1.0;
      for(int i=0; i<8; i++)
        {
         xnum = (xnum + p[i]) * z;
         xden = xden * z + q[i];
        }
      res=xnum/xden+1.0;
      if(y1<y)
        {
         //--- adjust result for case  0.0 < argument < 1.0
         res=res/y1;
        }
      else
      if(y1>y)
        {
         //--- adjust result for case  2.0 < argument < 12.0
         for(int i=0; i<n; i++)
           {
            res=res*y;
            y=y+1.0;
           }
        }
     }
   else
     {
      //--- evaluate for argument <=12.0,
      if(y<=xbig)
        {
         ysq=y*y;
         sum=c[6];
         for(int i=0; i<6; i++)
           {
            sum=sum/ysq+c[i];
           }
         sum = sum/y - y + logsqrt2pi;
         sum = sum + (y-0.5)*MathLog(y);
         res = MathExp(sum);
        }
      else
        {
         return(QPOSINF);
        }
     }
//--- final adjustments and return
   if(parity==1)
      res=-res;
   if(fact!=1.0)
      res=fact/res;
//--- return result     
   return(res);
  }
//+------------------------------------------------------------------+
//| MathGammaLog                                                     |
//+------------------------------------------------------------------+
//| Comments from original FORTRAN code:                             |
//| http://www.netlib.org/specfun/algama                             |
//|                                                                  |
//| This routine calculates the LOG(GAMMA) function for a positive   |
//| real argument X. Computation is based on an algorithm outlined   |
//| in references 1 and 2. The program uses rational functions that  |
//| theoretically approximate LOG(GAMMA) to at least 18 significant  |
//| decimal digits.  The approximation for X > 12 is from reference  |
//| 3, while approximations for X < 12.0 are similar to those in     |
//| reference 1, but are unpublished.  The accuracy achieved depends |
//| on the arithmetic system, the compiler, the intrinsic functions, |
//| and proper selection of the machine-dependent constants.         |
//|                                                                  |
//| References:                                                      |
//| 1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for |
//|   the Natural Logarithm of the Gamma Function,' Math. Comp. 21,  |
//|   1967, pp. 198-203.                                             |
//| 2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,|
//|   1969.                                                          |
//| 3) Hart, Et. Al., Computer Approximations, Wiley and sons, New   |
//|   York, 1968.                                                    |
//| Authors: W. J. Cody and L. Stoltz                                |
//|          Argonne National Laboratory                             |
//|                                                                  |
//| https://gcc.gnu.org/ml/fortran/2007-11/msg00061/gamma.diff       |
//| For negative arguments (where netlib would return +Inf)          |
//| we use the following identity:                                   |
//| lgamma(y) = log(pi/(|y*sin(pi*y)|)) - lgamma(-y)                 |
//+------------------------------------------------------------------+
double MathGammaLog(const double x)
  {
//--- mathematical constants
   const double pnt68=0.6796875e0;
   const double sqrtpi=0.9189385332046727417803297e0;
//--- machine dependent parameters
   const double xbig=2.55E305;
   const double xinf=1.79E308;
   const double eps=2.22E-16;
   const double frtbig=2.25E76;
//--- numerator and denominator coefficients for rational minimax
//--- approximation over (0.5,1.5)
   const double d1=-5.772156649015328605195174e-1;
   const double p1[8]=
     {
      4.945235359296727046734888e0,2.018112620856775083915565e2,
      2.290838373831346393026739e3,1.131967205903380828685045e4,
      2.855724635671635335736389e4,3.848496228443793359990269e4,
      2.637748787624195437963534e4,7.225813979700288197698961e3
     };
   const double q1[8]=
     {
      6.748212550303777196073036e1,1.113332393857199323513008e3,
      7.738757056935398733233834e3,2.763987074403340708898585e4,
      5.499310206226157329794414e4,6.161122180066002127833352e4,
      3.635127591501940507276287e4,8.785536302431013170870835e3
     };
//--- numerator and denominator coefficients for rational minimax
//--- approximation over (1.5,4.0)
   const double d2=4.227843350984671393993777e-1;
   const double p2[8]=
     {
      4.974607845568932035012064e0,5.424138599891070494101986e2,
      1.550693864978364947665077e4,1.847932904445632425417223e5,
      1.088204769468828767498470e6,3.338152967987029735917223e6,
      5.106661678927352456275255e6,3.074109054850539556250927e6
     };
   const double q2[8]=
     {
      1.830328399370592604055942e2,7.765049321445005871323047e3,
      1.331903827966074194402448e5,1.136705821321969608938755e6,
      5.267964117437946917577538e6,1.346701454311101692290052e7,
      1.782736530353274213975932e7,9.533095591844353613395747e6
     };
//--- numerator and denominator coefficients for rational minimax
//--- approximation over (4.0,12.0)
   const double d4=1.791759469228055000094023e0;
   const double p4[8]=
     {
      1.474502166059939948905062e4,2.426813369486704502836312e6,
      1.214755574045093227939592e8,2.663432449630976949898078e9,
      2.940378956634553899906876e10,1.702665737765398868392998e11,
      4.926125793377430887588120e11,5.606251856223951465078242e11
     };
   const double q4[8]=
     {
      2.690530175870899333379843e3,6.393885654300092398984238e5,
      4.135599930241388052042842e7,1.120872109616147941376570e9,
      1.488613728678813811542398e10,1.016803586272438228077304e11,
      3.417476345507377132798597e11,4.463158187419713286462081e11
     };
//--- coefficients for minimax approximation over (12, inf)
   const double c[7]=
     {
      -1.910444077728e-03,8.4171387781295e-04,
      -5.952379913043012e-04,7.93650793500350248e-04,
      -2.777777777777681622553e-03,8.333333333333333331554247e-02,
      5.7083835261e-03
     };

   double y=x;
   double res=0.0;
   double corr=0;
   double xm1,xm2,xm4,xden,xnum;
   double ysq=0;

   if((y>0 && y<=xbig))
     {
      if(y<=eps)
        {
         res=-MathLog(y);
        }
      else
      if(y<=1.5)
        {
         //--- eps < x <= 1.5
         if(y<pnt68)
           {
            corr= -MathLog(y);
            xm1 = y;
           }
         else
           {
            corr=0.0;
            xm1 =(y-0.5)-0.5;
           }
         if((y<=0.5 || y>=pnt68))
           {
            xden = 1.0;
            xnum = 0;
            for(int i=0; i<8; i++)
              {
               xnum = xnum*xm1 + p1[i];
               xden = xden*xm1 + q1[i];
              }
            res=corr+(xm1 *(d1+xm1*(xnum/xden)));
           }
         else
           {
            xm2=(y-0.5)-0.5;
            xden = 1.0;
            xnum = 0.0;
            for(int i=0; i<8; i++)
              {
               xnum=xnum*xm2+p2[i];
               xden=xden*xm2+q2[i];
              }
            res=corr+xm2 *(d2+xm2*(xnum/xden));
           }
        }
      else
      if(y<=4.0)
        {
         //--- 1.5 < x .<= 4.0
         xm2=y-2.0;
         xden = 1.0;
         xnum = 0.0;
         for(int i=0; i<8; i++)
           {
            xnum = xnum*xm2 + p2[i];
            xden = xden*xm2 + q2[i];
           }
         res=xm2 *(d2+xm2*(xnum/xden));
        }
      else
      if(y<=12.0)
        {
         //--- 4.0 < x <= 12.0
         xm4=y-4.0;
         xden = -1.0;
         xnum = 0.0;
         for(int i=0; i<8; i++)
           {
            xnum = xnum*xm4 + p4[i];
            xden = xden*xm4 + q4[i];
           }
         res=d4+xm4*(xnum/xden);
        }
      else
        {
         //--- evaluate for argument>=12.0,
         res=0.0;
         if(y<=frtbig)
           {
            res=c[6];
            ysq=y*y;
            for(int i=0; i<6; i++)
              {
               res=res/ysq+c[i];
              }
           }
         res=res/y;
         corr= MathLog(y);
         res = res + sqrtpi - 0.5*corr;
         res = res + y*(corr-1.0);
        }
     }
   else
     {
      //--- return for bad arguments
      //--- res=QPOSINF;
      //--- for negative arguments (where netlib would return +Inf)
      //--- we use the following identity:
      //--- lgamma(y) = log(pi/(|y*sin(pi*y)|)) - lgamma(-y)
      //--- https://gcc.gnu.org/ml/fortran/2007-11/msg00061/gamma.diff
      if(y<0 && MathFloor(y)!=y)
        {
         //--- for abs(y) very close to zero, we use a series expansion to
         //--- the first order in y to avoid overflow.
         if(y>-1.e-100)
            res=-2*MathLog(MathAbs(y))-MathGammaLog(-y);
         else
            res=MathLog(M_PI/MathAbs(y*MathSin(M_PI*y)))-MathGammaLog(-y);
        }
      else
        {
         //--- negative integer values
         res=QPOSINF;
        }
     }
//--- return result
   return(res);
  }
//+------------------------------------------------------------------+
//| MathBeta                                                         |
//+------------------------------------------------------------------+
double MathBeta(const double a,const double b)
  {
   return MathExp(MathBetaLog(a,b));
  }
//+------------------------------------------------------------------+
//| MathBetaLog                                                      |
//+------------------------------------------------------------------+
double MathBetaLog(const double a,const double b)
  {
   return(MathGammaLog(a)+MathGammaLog(b)-MathGammaLog(a+b));
  }
//+------------------------------------------------------------------+
//| Incomplete beta function                                         |
//+------------------------------------------------------------------+
//| The incomplete beta function ratio is the probability that a     |
//| random variable from a beta distribution having parameters p     |
//| and q will be less than or equal to x.                           |
//|                                                                  |
//| Input:                                                           |
//| x : upper limit of integration.  x must be in (0,1) inclusive.   |
//| p : first beta distribution parameter.  p must be>0.0.           |
//| q : second beta distribution parameter. q must be>0.0.           |
//|                                                                  |
//| Returns the value of the incomplete Beta function ratio.         |
//|                                                                  |
//| Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.     |
//| C version by John Burkardt.                                      |
//|                                                                  |
//| Reference:                                                       |
//|                                                                  |
//| KL Majumder, GP Bhattacharjee,                                   |
//| Algorithm AS 63: The incomplete Beta Integral,                   |
//| Applied Statistics,.Volume 22, Number 3, 1973, pages 409-411.    |
//+------------------------------------------------------------------+
double MathBetaIncomplete(const double x,const double p,const double q)
  {
   double acu=0.1E-16;
   double pp,qq,xx;
   int    ifault=0;
   double value=x;
//--- check the input arguments
   if(p<=0.0 || q<=0.0)
     {
      ifault=1;
      return(value);
     }
//--- check x range
   if(x<0.0 || 1.0<x)
     {
      ifault=2;
      return(value);
     }
//--- special cases
   if(x==0.0 || x==1.0)
     {
      return(value);
     }
//--- calc beta log
   double beta_log=MathGammaLog(p)+MathGammaLog(q)-MathGammaLog(p+q);
//--- change tail if necessary and determine S
   double psq= p+q;
   double cx = 1.0-x;

   int indx=0;
   if(p<psq*x)
     {
      xx = cx;
      cx = x;
      pp = q;
      qq = p;
      indx=1;
     }
   else
     {
      xx = x;
      pp = p;
      qq = q;
      indx=0;
     }

   value=1.0;
   double term=1.0;
   double ai=1.0;
   int    ns=(int)(qq+cx*psq);
//--- use the Soper reduction formula
   double rx=xx/cx;
   double temp=qq-ai;

   if(ns==0.0)
     {
      rx=xx;
     }

   for(;;)
     {
      term=term*temp*rx/(pp+ai);
      value= value+term;;
      temp = MathAbs(term);
      //---
      if(temp<=acu && temp<=acu*value)
        {
         value=value*MathExp(pp*MathLog(xx)+(qq-1.0)*MathLog(cx)-beta_log)/pp;
         if(indx)
           {
            value=1.0-value;
           }
         break;
        }

      ai=ai+1.0;
      ns--;

      if(0<=ns)
        {
         temp=qq-ai;
         if(ns==0)
           {
            rx=xx;
           }
        }
      else
        {
         temp= psq;
         psq = psq+1.0;
        }
     }
//--- return result
   return(value);
  }
//+------------------------------------------------------------------+
//| Incomplete Gamma function                                        |
//+------------------------------------------------------------------+
//| Inputs:                                                          |
//| x     : value at which the c.d.f is to be computed               |
//| alpha : parameter of the gamma function (>0)                     |
//+------------------------------------------------------------------+
//| Comment from original FORTRAN code:                              |
//| www.nist.gov/sites/default/files/documents/itl/sed/stspac.f      |
//| ibid .../itl/sed/SED_Note_86-2-2.pdf                             |
//|                                                                  |
//| CDFGAM  Written by Charles p. Reeve, Statistical Engineering     |
//|         Division, National Bureau of Standards, Gaithersburg.    |
//|                                                                  |
//| Computing the cumulative distribution function of the Gamma      |
//| distribution (also known as the Incomplete Gamma ratio) to a     |
//| specified accuracy (truncation error in the infinite series).    |
//| the algorithm, described in Ref. 2, is a modification of         |
//| the algorithm of Ref. 1. Three features have been added:         |
//|                                                                  |
//|  1) A precise method of meeting the truncation accuracy,         |
//|  2) Computation of the upper tail area by decrementing alpha     |
//|     when that method is more efficient, and                      |
//|  3) A constant uflo >= the underflow limit on the computer.      |
//|                                                                  |
//| References:                                                      |
//|                                                                  |
//| 1) Lau, Chi-Leung, 'A Simple Series for the Incomplete Gamma     |
//|    Integral', Algorithm AS 147, Applied Statistics, vol. 29,     |
//|    No. 1, 1980, pp. 113-114.                                     |
//|                                                                  |
//| 2) Reeve, Charles p., 'An Algorithm for Computing the Gamma      |
//|    C.D.F. to a Specified Accuracy', Statistical Engineering      |
//|    Division,Note 86-2, October 1986.                             |
//+------------------------------------------------------------------+
double MathGammaIncomplete(double x,double alpha)
  {
   double eps=10e-20;
   int    iflag=0;
   bool   ll;
   int    imax=5000;
   double uflo=1.0e-100;
   double cdfx=0.0;
   int    k=0;
//--- check for validity of arguments alpha and eps
   if(alpha<=uflo || eps<=uflo)
     {
      iflag=1;
      return(QNaN);
     }

   iflag=0;
//--- check for special case of x
   if(x<=0.0)
      return(QNaN);
//--- evaluate the gamma p.d.f. and check for underflow
   double dx=double(x);

   double pdfl=double(alpha-1.0)*MathLog(dx)-dx-MathGammaLog(alpha);

   if(pdfl<MathLog(uflo))
     {
      if(x>=alpha) cdfx=1.0;
     }
   else
     {
      double p=alpha;
      double u=MathExp(pdfl);
      //--- determine whether to increment or decrement alpha (a.k.a. p)
      ll=true;
      if(x>=p)
        {
         k=int(p);

         if(p<=double(k))
            k=k-1;

         double eta=p-double(k);
         double bl=double((eta-1.0)*MathLog(dx)-dx-MathGammaLog(eta));
         ll=bl>MathLog(eps);
        }
      //---
      double epsx=eps/x;
      //---
      if(ll==true)
        {
         //--- increment p
         for(int i=0; i<=imax; i++)
           {
            if(u<=epsx*(p-x))
              {
               return(cdfx);
              }
            u=x*u/p;
            cdfx=cdfx+u;
            p=p+1.0;
           }
         iflag=2;
        }
      else
        {
         //--- decrement p
         for(int j=1; j<=k; j++)
           {
            p=p-1.0;

            if(u<=epsx*(x-p))
               break;

            cdfx=cdfx+u;
            u=p*u/x;
           }
         cdfx=1.0-cdfx;
        }
     }
//--- return result
   return(cdfx);
  }
//+------------------------------------------------------------------+
//| Computes the binomial coefficient C(n,k)=n!/(k!*(n-k)!)          |
//|
//| The value is calculated in such a way as to avoid overflow and   |
//| roundoff.  The calculation is done in integer arithmetic.        |
//|                                                                  |
//| Parameters:                                                      |
//| Input, int N, K, are the values of N and K.                      |
//| Output: the number of combinations of N things taken K at a time.|
//|                                                                  |
//| Author: John Burkardt                                            |
//|                                                                  |
//| Reference:                                                       |
//| ML Wolfson, HV Wright,                                           |
//| Algorithm 160: Combinatorial of M Things Taken N at a Time,      |
//| Communications of the ACM, Volume 6, N4, April 1963, page 161.   |
//+------------------------------------------------------------------+
long MathBinomialCoefficient(const int n,const int k)
  {
   int  mn,mx;
   long value=0;
   int n_k=n-k;
//---
   if(k<n_k)
     {
      mn = k;
      mx = n_k;
     }
   else
     {
      mn = n_k;
      mx = k;
     }
   if(mn>0)
     {
      value=mx+1;
      for(int i=2; i<=mn; i++)
         value=(value*(mx+i))/i;
     }
   else
   if(mn<0)
      return(0);
   else
   if(mn==0)
      return(1);
//---
   return(value);
  }
//+------------------------------------------------------------------+
//| MathBinomialCoefficientLog                                       |
//+------------------------------------------------------------------+
double MathBinomialCoefficientLog(const int n,const int k)
  {
   int  mn,mx;
   long value=0;
   int n_k=n-k;
//---
   if(k<n_k)
     {
      mn = k;
      mx = n_k;
     }
   else
     {
      mn = n_k;
      mx = k;
     }
   if(mn>0)
     {
      value=mx+1;
      for(int i=2; i<=mn; i++)
         value=(value*(mx+i))/i;
     }
   else
   if(mn<0)
      return(QNEGINF);
   else
   if(mn==0)
      return(0);
//---
   return MathLog(value);
  }
//+------------------------------------------------------------------+
//| Computes the logarithm of the Binomial coefficient.              |
//|                                                                  |
//| Log(C(n,k))=Log(n!/(k!*(n-k)!))                                  |
//|                                                                  |
//| Parameters:                                                      |
//| Input, int n, k, are the values of n and k.                      |
//| Return value: the logarithm of C(n,k).                           |
//+------------------------------------------------------------------+
double MathBinomialCoefficientLog(const double n,const double k)
  {
   return(MathGammaLog(n+1)-MathGammaLog(k+1)-MathGammaLog(n-k+1));
  }
//+------------------------------------------------------------------+
//| Computes the function Hypergeometric_2F2(a,b,c,d,z) using        |
//| a Taylor method.                                                 |
//|                                                                  |
//| Author    : John Pearson                                         |
//| Reference : MSc thesis "Computation of Hypergeometric Functions" |
//+------------------------------------------------------------------+
double MathHypergeometric2F2(const double a,const double b,const double c,const double d,const double z)
  {
   double a1=1.0;
   double b1=1.0;
   double tol=10E-10;  // set tolerance
//--- direct summation
   for(int j=1; j<=500; j++)
     {
      //--- update current term in terms of previous one
      a1=(a+j-1)*(b+j-1)/(c+j-1)/(d+j-1)*z/j*a1;
      //--- update sum of terms computed so far
      b1=b1+a1;
      //--- terminate summation if stopping criterion is satisfied
      if(MathAbs(a1)/MathAbs(b1)<tol)
         return(b1);
     }
//--- unknown
   return(QNaN);
  }
//+------------------------------------------------------------------+
//| Returns log(0) depending on tail and log_mode arugments          |
//+------------------------------------------------------------------+
double TailLog0(const bool tail,const bool log_mode)
  {
   if(tail==true)
     {
      if(log_mode==true)
         return(QNEGINF);
      else
         return(0.0);
     }
   else
     {
      if(log_mode==true)
         return(0.0);
      else
         return(1.0);
     }
  }
//+------------------------------------------------------------------+
//| Returns log(1) depending on tail and log_mode arugments          |
//+------------------------------------------------------------------+
double TailLog1(const bool tail,const bool log_mode)
  {
   if(tail==true)
     {
      if(log_mode==true)
         return(0.0);
      else
         return(1.0);
     }
   else
     {
      if(log_mode==true)
         return(1.0);
      else
         return(0.0);
     }
  }
//+------------------------------------------------------------------+
//| Returns value depending on tail and log_mode arugments           |
//+------------------------------------------------------------------+
double TailLogValue(const double value,const bool tail,const bool log_mode)
  {
   if(tail==true)
     {
      if(log_mode)
         return MathLog(MathAbs(value));
      else
         return(value);
     }
   else
     {
      if(log_mode)
         return MathLog(MathAbs(1.0-value));
      else
         return(1.0-value);
     }
  }
//+------------------------------------------------------------------+
//| Returns value depending on tail and log_mode arugments           |
//+------------------------------------------------------------------+
double TailLogProbability(const double probability,const bool tail,const bool log_mode)
  {
   if(tail==true)
     {
      if(log_mode)
         return MathExp(probability);
      else
         return(probability);
     }
   else
     {
      if(log_mode)
         return(1.0-MathExp(probability));
      else
         return(1.0-probability);
     }
  }
//+------------------------------------------------------------------+
//| MathSequence                                                     |
//+------------------------------------------------------------------+
//| The function generates a sequence of double numbers              |
//| with given starting and ending values and step.                  |
//|                                                                  |
//| Arguments:                                                       |
//| from        : Starting value of the sequence                     |
//| to          : Last value of the sequence                         |
//| step        : Step of the sequence                               |
//| result[]    : Array for the sequence values                      |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSequence(const double from,const double to,const double step,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(from) || !MathIsValidNumber(to) || !MathIsValidNumber(step))
      return(false);

   if(to<from)
      return(false);
//--- calculate number of elements and prepare output array
   int count=1+int((to-from)/step);

   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare sequence
   for(int i=0; i<count; i++)
      result[i]=from+i*step;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSequence                                                     |
//+------------------------------------------------------------------+
//| The function generates a sequence of integer numbers             |
//| with given starting and ending values and step.                  |
//|                                                                  |
//| Arguments:                                                       |
//| from        : Starting value of the sequence                     |
//| to          : Last value of the sequence                         |
//| step        : Step of the sequence                               |
//| result[]    : Array for the sequence values                      |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSequence(const int from,const int to,const int step,int &result[])
  {
   if(to<from)
      return(false);
//--- calculate number of elements and prepare output array
   int count=1+int((to-from)/step);

   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare sequence
   for(int i=0; i<count; i++)
      result[i]=from+i*step;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSequenceByCount                                              |
//+------------------------------------------------------------------+
//| The function generates a sequence of double numbers              |
//| with given starting and ending values and number of elements     |
//| in the sequence.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| from        : Starting value of the sequence                     |
//| to          : Last value of the sequence                         |
//| count       : Number of elements in the sequence                 |
//| result[]    : Array for sequence values                          |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSequenceByCount(const double from,const double to,const int count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(from) || !MathIsValidNumber(to))
      return(false);

   if(to<from || count<=0)
      return(false);
//--- prepare output array
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare sequence
   double step=(to-from)/(count-1);

   for(int i=0; i<count; i++)
      result[i]=from+i*step;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSequenceByCount                                              |
//+------------------------------------------------------------------+
//| The function generates a sequence of double numbers              |
//| with given starting and ending values and number of elements     |
//| in the sequence.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| from        : Starting value of the sequence                     |
//| to          : Last value of the sequence                         |
//| count       : Number of elements in the sequence                 |
//| result[]    : Array for sequence values                          |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSequenceByCount(const int from,const int to,const int count,int &result[])
  {
   if(to<from || count<=0)
      return(false);
//--- prepare output array
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare sequence
   int step=(to-from)/(count-1);
   for(int i=0; i<count; i++)
      result[i]=from+i*step;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReplicate                                                    |
//+------------------------------------------------------------------+
//| The function replicates the sequence of double numbers           |
//| from array[] with given times (count).                           |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values to replicate                     |
//| count       : Number of elements in the sequence                 |
//| result[]    : Array for sequence values                          |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReplicate(const double &array[],const int count,double &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- check replicate count
   if(count<=0)
      return(false);

//--- prepare output array
   int target_size=size*count;
   if(ArraySize(result)<target_size)
      if(ArrayResize(result,target_size)!=target_size)
         return(false);

//--- replicate array elements
   for(int i=0; i<count; i++)
      ArrayCopy(result,array,i*size,0,WHOLE_ARRAY);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReplicate                                                    |
//+------------------------------------------------------------------+
//| The function replicates the sequence of integer numbers          |
//| from array[] with given times (count).                           |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values to replicate                      |
//| count       : Number of elements in the sequence                 |
//| result[]    : Array for sequence values                          |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReplicate(const int &array[],const int count,int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- check replicate count
   if(count<=0)
      return(false);

//--- prepare output array
   int target_size=size*count;
   if(ArraySize(result)<target_size)
      if(ArrayResize(result,target_size)!=target_size)
         return(false);

//--- replicate array elements
   for(int i=0; i<count; i++)
      ArrayCopy(result,array,i*size,0,WHOLE_ARRAY);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReverse                                                      |
//+------------------------------------------------------------------+
//| The function reverses order of the elements from the array[].    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for reversed elements                        |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReverse(const double &array[],double &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate elements of the reversed array
   for(int i=0; i<size; i++)
     {
      int idx=size-1-i;
      result[idx]=array[i];
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReverse                                                      |
//+------------------------------------------------------------------+
//| The function reverses order of the elements from the array[].    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for reversed elements                        |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReverse(const int &array[],int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate elements of the reversed array
   for(int i=0; i<size; i++)
     {
      int idx=size-1-i;
      result[idx]=array[i];
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReverse                                                      |
//+------------------------------------------------------------------+
//| The function reverses order of the elements from the array[].    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReverse(double &array[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
   int count=size/2;
   for(int i=0; i<count; i++)
     {
      int idx=size-1-i;
      //--- swap array[i] and array[size-1-i];
      double t=array[i];
      array[i]=array[idx];
      array[idx]=t;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathReverse                                                      |
//+------------------------------------------------------------------+
//| The function reverses order of the elements from the array[].    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathReverse(int &array[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
   int count=size/2;
   for(int i=0; i<count; i++)
     {
      int idx=size-1-i;
      //--- swap array[i] and array[size-1-i];
      int t=array[i];
      array[i]=array[idx];
      array[idx]=t;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathIdentical                                                    |
//+------------------------------------------------------------------+
//| The function compares two arrays.                                |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array                                        |
//| array2[]    : Second array                                       |
//|                                                                  |
//| Return value: true if the arrays are equal, otherwise false.     |
//+------------------------------------------------------------------+
bool MathIdentical(const double &array1[],const double &array2[])
  {
   int size1=ArraySize(array1);
   int size2=ArraySize(array1);
//--- check size of the arrays
   if(size1!=size2)
      return(false);
//--- check empty array case
   if(size1==0)
      return(false);
//--- check all array elements
   for(int i=0; i<size1; i++)
     {
      if(array1[i]!=array2[i])
         return(false);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathIdentical                                                    |
//+------------------------------------------------------------------+
//| The function compares two arrays.                                |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array                                        |
//| array2[]    : Second array                                       |
//|                                                                  |
//| Return value: true if the arrays are equal, otherwise false.     |
//+------------------------------------------------------------------+
bool MathIdentical(const int &array1[],const int &array2[])
  {
   int size1=ArraySize(array1);
   int size2=ArraySize(array1);
//--- check size of the arrays
   if(size1!=size2)
      return(false);
//--- check empty array case
   if(size1==0)
      return(false);
//--- check all array elements
   for(int i=0; i<size1; i++)
     {
      if(array1[i]!=array2[i])
         return(false);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathUnique                                                       |
//+------------------------------------------------------------------+
//| The function extracts the unique values from the array.          |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for unique values                            |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathUnique(const double &array[],double &result[])
  {
//--- check array size
   int size=ArraySize(array);
   if(size==0)
      return(false);

//--- prepare additional array
   bool checked[];
   if(ArrayResize(checked,size)!=size)
      return(false);
   ArrayFill(checked,0,size,false);

//--- prepare result array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- find unique elements
   int unique_count=0;
   double value=0;
   for(;;)
     {
      bool flag=false;
      for(int i=unique_count; i<size; i++)
        {
         if(!flag && !checked[i])
           {
            value=array[i];
            result[unique_count]=array[i];
            unique_count++;
            checked[i]=true;
            flag=true;
           }
         else
         if(flag && value==array[i])
            checked[i]=true;
        }
      if(!flag)
         break;
     }
//--- resize target array
   ArrayResize(result,unique_count);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathUnique                                                       |
//+------------------------------------------------------------------+
//| The function extracts the unique values from the array.          |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for unique values                            |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathUnique(const int &array[],int &result[])
  {
//--- check array size
   int size=ArraySize(array);
   if(size==0)
      return(false);

//--- prepare additional array
   bool checked[];
   if(ArrayResize(checked,size)!=size)
      return(false);
   ArrayFill(checked,0,size,false);

//--- prepare result array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- find unique elements
   int unique_count=0;
   int value=0;
   for(;;)
     {
      bool flag=false;
      for(int i=unique_count; i<size; i++)
        {
         if(!flag && !checked[i])
           {
            value=array[i];
            result[unique_count]=array[i];
            unique_count++;
            checked[i]=true;
            flag=true;
           }
         else
         if(flag && value==array[i])
            checked[i]=true;
        }
      if(!flag)
         break;
     }
//--- resize target array
   ArrayResize(result,unique_count);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathQuickSortAscending                                           |
//+------------------------------------------------------------------+
//| The function sorts array[] and indices[] simultaneously          |
//| using QuickSort algorithm. After sorting the initial ordering    |
//| of the elements is located in the indices[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| array       : Array with values to sort                          |
//| indices     : Targe array to sort                                |
//| first       : First element index                                |
//| last        : Last element index                                 |
//|                                                                  |
//| Return value: None                                               |
//+------------------------------------------------------------------+
void MathQuickSortAscending(double &array[],int &indices[],int first,int last)
  {
   int    i,j,t_int;
   double p_double,t_double;
//--- check
   if(first<0 || last<0)
      return;
//--- sort
   i=first;
   j=last;
   while(i<last)
     {
      //--- ">>1" is quick division by 2
      p_double=array[(first+last)>>1];
      while(i<j)
        {
         while(array[i]<p_double)
           {
            //--- control the output of the array bounds
            if(i==ArraySize(array)-1)
               break;
            i++;
           }
         while(array[j]>p_double)
           {
            //--- control the output of the array bounds
            if(j==0)
               break;
            j--;
           }
         if(i<=j)
           {
            //-- swap elements i and j
            t_double=array[i];
            array[i]=array[j];
            array[j]=t_double;
            //-- swap indices i and j
            t_int=indices[i];
            indices[i]=indices[j];
            indices[j]=t_int;
            i++;
            //--- control the output of the array bounds
            if(j==0)
               break;
            j--;
           }
        }
      if(first<j)
         MathQuickSortAscending(array,indices,first,j);
      first=i;
      j=last;
     }
  }
//+------------------------------------------------------------------+
//| MathQuickSortDescending                                          |
//+------------------------------------------------------------------+
//| The function sorts array[] and indices[] simultaneously          |
//| using QuickSort algorithm. After sorting the initial ordering    |
//| of the elements is located in the indices[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| array       : Array with values to sort                          |
//| indices     : Targe array to sort                                |
//| first       : First element index                                |
//| last        : Last element index                                 |
//|                                                                  |
//| Return value: None                                               |
//+------------------------------------------------------------------+
void MathQuickSortDescending(double &array[],int &indices[],int first,int last)
  {
   int    i,j,t_int;
   double p_double,t_double;
//--- check
   if(first<0 || last<0)
      return;
//--- sort
   i=first;
   j=last;
   while(i<last)
     {
      //--- ">>1" is quick division by 2
      p_double=array[(first+last)>>1];
      while(i<j)
        {
         while(array[i]>p_double)
           {
            //--- control the output of the array bounds
            if(i==ArraySize(array)-1)
               break;
            i++;
           }
         while(array[j]<p_double)
           {
            //--- control the output of the array bounds
            if(j==0)
               break;
            j--;
           }
         if(i<=j)
           {
            //-- swap elements i and j
            t_double=array[i];
            array[i]=array[j];
            array[j]=t_double;
            //-- swap indices i and j
            t_int=indices[i];
            indices[i]=indices[j];
            indices[j]=t_int;
            i++;
            //--- control the output of the array bounds
            if(j==0)
               break;
            j--;
           }
        }
      if(first<j)
         MathQuickSortDescending(array,indices,first,j);
      first=i;
      j=last;
     }
  }
//+------------------------------------------------------------------+
//| MathQuickSort                                                    |
//+------------------------------------------------------------------+
//| The function sorts array[] and indices[] simultaneously          |
//| using QuickSort algorithm. After sorting the initial ordering    |
//| of the elements is located in the indices[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| array       : Array with values to sort                          |
//| indices     : Targe array to sort                                |
//| first       : First element index                                |
//| last        : Last element index                                 |
//| mode        : Sort direction (>0 ascending,otherwise descending) |
//|                                                                  |
//| Return value: None                                               |
//+------------------------------------------------------------------+
void MathQuickSort(double &array[],int &indices[],int first,int last,int mode)
  {
   if(mode>0)
      MathQuickSortAscending(array,indices,first,last);
   else
      MathQuickSortDescending(array,indices,first,last);
  }
//+------------------------------------------------------------------+
//| MathOrder                                                        |
//+------------------------------------------------------------------+
//| The function calculates the order of elements from the array.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for sorted indices                           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathOrder(const double &array[],int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- prepare temporary array
   double tmp_array[];
   ArrayCopy(tmp_array,array,0,0,WHOLE_ARRAY);
//--- prepare identical permutation
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
      result[i]=i+1;
//--- calculate order of numbers
   MathQuickSortAscending(tmp_array,result,0,size-1);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathOrder                                                        |
//+------------------------------------------------------------------+
//| The function calculates the order of elements from the array.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Array for sorted indices                           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathOrder(const int &array[],int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return(false);
//--- prepare temporary array
   double tmp_array[];
   if(ArrayResize(tmp_array,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      tmp_array[i]=array[i];
//--- prepare identical permutation
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
      result[i]=i+1;
//--- calculate order of numbers
   MathQuickSortAscending(tmp_array,result,0,size-1);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseNot                                                   |
//+------------------------------------------------------------------+
//| The function calculates NOT logical operation for elements       |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseNot(const int &array[],int &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);

//--- prepare target array and calculate ~array[i]
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
      result[i]=~array[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseNot                                                   |
//+------------------------------------------------------------------+
//| The function calculates NOT logical operation for elements       |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseNot(int &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate ~array[i]
   for(int i=0; i<size; i++)
      array[i]=~array[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseAnd                                                   |
//+------------------------------------------------------------------+
//| The function calculates AND logical operation for elements       |
//| of the array1[] and array2[].                                    |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with values                            |
//| array2[]    : Second array with values                           |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseAnd(const int &array1[],const int &array2[],int &result[])
  {
   int size1=ArraySize(array1);
   int size2=ArraySize(array2);
   if(size1==0 || size2==0 || size1!=size2)
      return(false);
//--- prepare target array and calculate array1&array2
   if(ArraySize(result)<size1)
      if(ArrayResize(result,size1)!=size1)
         return(false);

   for(int i=0; i<size1; i++)
      result[i]=array1[i]&array2[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseOr                                                    |
//+------------------------------------------------------------------+
//| The function calculates OR logical operation for elements        |
//| of the array1[] and array2[].                                    |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with values                            |
//| array2[]    : Second array with values                           |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseOr(const int &array1[],const int &array2[],int &result[])
  {
   int size1=ArraySize(array1);
   int size2=ArraySize(array2);
   if(size1==0 || size2==0 || size1!=size2)
      return(false);
//--- prepare target array and calculate array1|array2
   if(ArraySize(result)<size1)
      if(ArrayResize(result,size1)!=size1)
         return(false);

   for(int i=0; i<size1; i++)
      result[i]=array1[i]|array2[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseXor                                                   |
//+------------------------------------------------------------------+
//| The function calculates XOR logical operation for elements       |
//| of the array1[] and array2[].                                    |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with values                            |
//| array2[]    : Second array with values                           |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseXor(const int &array1[],const int &array2[],int &result[])
  {
   int size1=ArraySize(array1);
   int size2=ArraySize(array2);
   if(size1==0 || size2==0 || size1!=size2)
      return(false);
//--- prepare target array and calculate array1^array2
   if(ArraySize(result)<size1)
      if(ArrayResize(result,size1)!=size1)
         return(false);

   for(int i=0; i<size1; i++)
      result[i]=array1[i]^array2[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseShiftL                                                |
//+------------------------------------------------------------------+
//| The function calculates << logical operation for elements        |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| n           : Number of bits to shift                            |
//| result[]    : Array for sorted indices                           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseShiftL(const int &array[],const int n,int &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate array<<n
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   for(int i=0; i<size; i++)
      result[i]=array[i]<<n;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseShiftL                                                |
//+------------------------------------------------------------------+
//| The function calculates << logical operation for elements        |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| n           : Number of bits to shift                            |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseShiftL(int &array[],const int n)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate array<<n
   for(int i=0; i<size; i++)
      array[i]=array[i]<<n;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseShiftR                                                |
//+------------------------------------------------------------------+
//| The function calculates >> logical operation for elements        |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| n           : Number of bits to shift                            |
//| result[]    : Array for sorted indices                           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseShiftR(const int &array[],const int n,int &result[])
  {
   int size=ArraySize(array);
//---
   if(size==0)
      return(false);
//--- prepare target array and calculate array>>n
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   for(int i=0; i<size; i++)
      result[i]=array[i]>>n;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathBitwiseShiftR                                                |
//+------------------------------------------------------------------+
//| The function calculates >> logical operation for elements        |
//| of the array[].                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| n           : Number of bits to shift                            |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathBitwiseShiftR(int &array[],const int n)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate array>>n
   for(int i=0; i<size; i++)
      array[i]=array[i]>>n;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeSum                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative sum of the elements           |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with cumulative sum                   |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeSum(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate cumulative sum
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   double sum=0.0;
   for(int i=0; i<size; i++)
     {
      //--- check value
      if(!MathIsValidNumber(array[i]))
         return(false);
      sum+=array[i];
      result[i]=sum;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeSum                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative sum of the elements           |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeSum(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate cumulative sum
   double sum=0.0;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]))
         return(false);
      sum+=array[i];
      array[i]=sum;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeProduct                                            |
//+------------------------------------------------------------------+
//| The function calculates cumulative product of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with cumulative product               |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeProduct(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate cumulative product
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   double product=1.0;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]))
         return(false);
      product*=array[i];
      result[i]=product;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeProduct                                            |
//+------------------------------------------------------------------+
//| The function calculates cumulative product of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeProduct(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate cumulative product
   double product=1.0;
   for(int i=0; i<size; i++)
     {
      //--- check value
      if(!MathIsValidNumber(array[i]))
         return(false);
      product*=array[i];
      array[i]=product;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeMin                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative minimum of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with cumulative product               |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeMin(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate cumulative min
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   double min=QPOSINF;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]) || !MathIsValidNumber(min))
         min+=array[i];
      else
         min=MathMin(min,array[i]);
      result[i]=min;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeMin                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative minimum of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeMin(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);

   double min=QPOSINF;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]) || !MathIsValidNumber(min))
         min+=array[i];
      else
         min=MathMin(min,array[i]);
      array[i]=min;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeMax                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative maximum of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with cumulative product               |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeMax(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate cumulative max
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);

   double max=QNEGINF;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]) || !MathIsValidNumber(max))
         max+=array[i];
      else
         max=MathMax(max,array[i]);
      result[i]=max;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeMax                                                |
//+------------------------------------------------------------------+
//| The function calculates cumulative maximum of the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeMax(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);

   double max=QNEGINF;
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]) || !MathIsValidNumber(max))
         max+=array[i];
      else
         max=MathMax(max,array[i]);
      array[i]=max;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSin                                                          |
//+------------------------------------------------------------------+
//| The function calculates sin(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSin(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathSin(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSin                                                          |
//+------------------------------------------------------------------+
//| The function calculates sin(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSin(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathSin(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCos                                                          |
//+------------------------------------------------------------------+
//| The function calculates cos(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCos(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathCos(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCos                                                          |
//+------------------------------------------------------------------+
//| The function calculates cos(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCos(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathCos(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTan                                                          |
//+------------------------------------------------------------------+
//| The function calculates tan(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTan(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathTan(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTan                                                          |
//+------------------------------------------------------------------+
//| The function calculates tan(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTan(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathTan(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArcsin                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arcsin(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArcsin(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArcsin(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArcsin                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arcsin(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArcsin(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArcsin(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArccos                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arccos(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArccos(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArccos(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArccos                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arccos(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArccos(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArccos(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArctan                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arctan(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArctan(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArctan(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArctan                                                       |
//+------------------------------------------------------------------+
//| The function calculates Arctan(x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArctan(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArctan(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSinPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates sin(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSinPi(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathSin(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSinPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates sin(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSinPi(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathSin(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCosPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates cos(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCosPi(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathCos(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCosPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates cos(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCosPi(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathCos(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTanPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates tan(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTanPi(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathTan(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTanPi                                                        |
//+------------------------------------------------------------------+
//| The function calculates tan(pi*x) for the elements               |
//| from the array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTanPi(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathTan(M_PI*array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathAbs                                                          |
//+------------------------------------------------------------------+
//| The function calculates abs(x) for the elements from the array.  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathAbs(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathAbs(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathAbs                                                          |
//+------------------------------------------------------------------+
//| The function calculates abs(x) for the elements from the array.  |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathAbs(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathAbs(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCeil                                                         |
//+------------------------------------------------------------------+
//| The function calculates integer numberic values closest from     |
//| above for the elements from the array[].                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCeil(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate Ceil(array[i])
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
      result[i]=MathCeil(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCeil                                                         |
//+------------------------------------------------------------------+
//| The function calculates integer numberic values closest from     |
//| above for the elements from the array[].                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCeil(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
   for(int i=0; i<size; i++)
      array[i]=MathCeil(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathFloor                                                        |
//+------------------------------------------------------------------+
//| The function calculates integer numberic values closest from     |
//| below for the elements from the array[].                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathFloor(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate Floor(array[i])
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
      result[i]=MathFloor(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathFloor                                                        |
//+------------------------------------------------------------------+
//| The function calculates integer numberic values closest from     |
//| below for the elements from the array[].                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathFloor(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
   for(int i=0; i<size; i++)
      array[i]=MathFloor(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTrunc                                                        |
//+------------------------------------------------------------------+
//| The function calculates trunc(x) for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTrunc(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array and calculate Trunc(array[i])
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
   for(int i=0; i<size; i++)
     {
      if(array[i]>=0)
         result[i]=MathFloor(array[i]);
      else
         result[i]=MathCeil(array[i]);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTrunc                                                        |
//+------------------------------------------------------------------+
//| The function calculates trunc(x) for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTrunc(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
   for(int i=0; i<size; i++)
     {
      if(array[i]>=0)
         array[i]=MathFloor(array[i]);
      else
         array[i]=MathCeil(array[i]);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSqrt                                                         |
//+------------------------------------------------------------------+
//| The function calculates sqrt(x) for the elements from the array. |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSqrt(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathSqrt(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSqrt                                                         |
//+------------------------------------------------------------------+
//| The function calculates sqrt(x) for the elements from the array. |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSqrt(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathSqrt(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathExp                                                          |
//+------------------------------------------------------------------+
//| The function calculates exp(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathExp(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathExp(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathExp                                                          |
//+------------------------------------------------------------------+
//| The function calculates exp(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathExp(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathExp(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathPow                                                          |
//+------------------------------------------------------------------+
//| The function calculates Pow(x,power) for the elements            |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathPow(const double &array[],const double power,double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NAN
   if(!MathIsValidNumber(power))
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathPow(array[i],power);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathPow                                                          |
//+------------------------------------------------------------------+
//| The function calculates Pow(x,power) for the elements            |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathPow(double &array[],const double power)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NAN
   if(!MathIsValidNumber(power))
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathPow(array[i],power);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog                                                          |
//+------------------------------------------------------------------+
//| The function calculates log(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathLog(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog                                                          |
//+------------------------------------------------------------------+
//| The function calculates log(x) for the elements from the array[].|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathLog(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog                                                          |
//+------------------------------------------------------------------+
//| The function calculates logarithm with the specified base        |
//| for the elements from the array[].                               |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog(const double &array[],const double base,double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NAN
   if(!MathIsValidNumber(base))
      return(false);
//--- check base
   if(base==1.0 || base<=0.0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   double factor=1.0/MathLog(base);
   for(int i=0; i<size; i++)
      result[i]=MathLog(array[i])*factor;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog                                                          |
//+------------------------------------------------------------------+
//| The function calculates logarithm with the specified base        |
//| for the elements from the array[].                               |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog(double &array[],const double base)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NAN
   if(!MathIsValidNumber(base))
      return(false);
//--- check base
   if(base==1.0 || base<=0.0)
      return(false);
//--- calculate values
   double factor=1.0/MathLog(base);
   for(int i=0; i<size; i++)
      array[i]=MathLog(array[i])*factor;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog2                                                         |
//+------------------------------------------------------------------+
//| The function calculates logarithms on base 2 for the elements    |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog2(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   static const double factor=1.0/M_LN2;
   for(int i=0; i<size; i++)
      result[i]=MathLog(array[i])*factor;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog2                                                         |
//+------------------------------------------------------------------+
//| The function calculates logarithms on base 2 for the elements    |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog2(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   static const double factor=1.0/M_LN2;
   for(int i=0; i<size; i++)
      array[i]=MathLog(array[i])*factor;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog10                                                        |
//+------------------------------------------------------------------+
//| The function calculates logarithm on base 10 for the elements    |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog10(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathLog10(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog10                                                        |
//+------------------------------------------------------------------+
//| The function calculates logarithm on base 10 for the elements    |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog10(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathLog10(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArctan2                                                      |
//+------------------------------------------------------------------+
//| The function calculates polar angles of vector (x[i],y[i])       |
//| in range [-pi, pi] for the elements from the x[] and y[] arrays. |
//|                                                                  |
//| Arguments:                                                       |
//| x[]         : Array with double values                           |
//| y[]         : Array with double values                           |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArctan2(const double &x[],const double &y[],double &result[])
  {
   int size=ArraySize(x);
   if(size==0 || ArraySize(y)!=size)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArctan2(x[i],y[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathRound                                                        |
//+------------------------------------------------------------------+
//| The function calculates rounded values to a given precision      |
//| for the elements from the array[].                               |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| digits      : Precision                                          |
//| result[]    : Output array with rounded values                   |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathRound(const double &array[],int digits,double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
   if(digits<0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathRound(array[i],digits);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathRound                                                        |
//+------------------------------------------------------------------+
//| The function calculates rounded values to a given precision      |
//| for the elements from the array[].                               |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| digits      : Precision                                          |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathRound(double &array[],int digits)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
   if(digits<0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathRound(array[i],digits);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathDifference                                                   |
//+------------------------------------------------------------------+
//| The function calculates laggged differences of the elements      |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| lag         : Lag value                                          |
//| result[]    : Output array with calculated differences           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathDifference(const double &array[],const int lag,double &result[])
  {
   int size=ArraySize(array)-lag;
   if(lag<1 || size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=array[i+lag]-array[i];
//--- change array size for call with differences
   if(ArrayResize(result,size)!=size)
      return(false);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathDifference                                                   |
//+------------------------------------------------------------------+
//| The function calculates laggged differences of the elements      |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with integer values                          |
//| lag         : Lag value                                          |
//| result[]    : Output array with calculated differences           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathDifference(const int &array[],const int lag,int &result[])
  {
   int size=ArraySize(array)-lag;
   if(lag<1 || size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=array[i+lag]-array[i];
//--- change array size for call with differences
   if(ArrayResize(result,size)!=size)
      return(false);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathDifference                                                   |
//+------------------------------------------------------------------+
//| The function calculates laggged and iterated differences         |
//| of the elements from the array[].                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| lag         : Lag value                                          |
//| differences : Number of iterations                               |
//| result[]    : Output array with calculated differences           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathDifference(const double &array[],const int lag,const int differences,double &result[])
  {
   int size=ArraySize(array);
//--- check parameters
   if(size==0 || differences<0 || lag*differences>=size)
      return(false);
//--- copy data from initial array
   if(ArrayCopy(result,array,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate iterated differences
   for(int i=0; i<differences; i++)
     {
      if(!MathDifference(result,lag,result))
         return false;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathDifference                                                   |
//+------------------------------------------------------------------+
//| The function calculates laggged and iterated differences         |
//| of the elements from the array[].                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with integer values                          |
//| lag         : Lag value                                          |
//| differences : Number of iterations                               |
//| result[]    : Output array with calculated differences           |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathDifference(const int &array[],const int lag,const int differences,int &result[])
  {
   int size=ArraySize(array);
//--- check parameters
   if(size==0 || differences<0 || lag*differences>=size)
      return(false);
//--- copy data from initial array
   if(ArrayCopy(result,array,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate iterated differences
   for(int i=0; i<differences; i++)
     {
      if(!MathDifference(result,lag,result))
         return false;
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[].                                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| count       : Number of elements to choose                       |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const double &array[],const int count,double &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return false;
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
   for(int i=0; i<count; i++)
     {
      int ind=(int)(count-1)*MathRand()/32767;
      result[i]=array[ind];
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[].                                         |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with integer values                          |
//| count       : Number of elements to choose                       |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const int &array[],const int count,int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0)
      return false;
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
   for(int i=0; i<count; i++)
     {
      int ind=(int)(count-1)*MathRand()/32767;
      result[i]=array[ind];
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] with replacement.                        |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| count       : Number of elements to choose                       |
//| replace     : If true, only unique values are placed to result[] |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const double &array[],const int count,const bool replace,double &result[])
  {
//--- unique values not needed
   if(!replace)
      return MathSample(array,count,result);
   int size=ArraySize(array);
   if(size==0 || count>size)
      return(false);
//--- prepare target array
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- unique values needed, prepare indices
   int indices[];
   MathSequenceByCount(0,count-1,count,indices);
   int swaps=size/2;
   int max_v=count-1;
   for(int i=0; i<swaps; i++)
     {
      //--- select random indices and swap them
      int ind1=(int)max_v*MathRand()/32767;
      int ind2=(int)max_v*MathRand()/32767;
      int t=indices[ind1];
      indices[ind1]=indices[ind2];
      indices[ind2]=t;
     }
//--- select data according to indices
   for(int i=0; i<count; i++)
      result[i]=array[indices[i]];
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] with replacement.                        |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with integer values                          |
//| count       : Number of elements to choose                       |
//| replace     : If true, only unique values are placed to result[] |
//| result[]    : Output array                                       |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const int &array[],const int count,const bool replace,int &result[])
  {
//--- unique values not needed
   if(!replace)
      return MathSample(array,count,result);
   int size=ArraySize(array);
   if(size==0 || count>size)
      return(false);
//--- prepare target array
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- unique values needed, prepare indices
   int indices[];
   MathSequenceByCount(0,count-1,count,indices);
   int swaps=size/2;
   int max_v=count-1;
   for(int i=0; i<swaps; i++)
     {
      //--- select random indices and swap them
      int ind1=(int)max_v*MathRand()/32767;
      int ind2=(int)max_v*MathRand()/32767;
      int t=indices[ind1];
      indices[ind1]=indices[ind2];
      indices[ind2]=t;
     }
//--- select data according to indices
   for(int i=0; i<count; i++)
      result[i]=array[indices[i]];
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] according to the given probabilities.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]         : Array with double values                       |
//| probabilities[] : Array with probabilities                       |
//| count           : Number of elements to choose                   |
//| result[]        : Output array                                   |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const double &array[],double &probabilities[],const int count,double &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0 || size!=ArraySize(probabilities))
      return(false);

   if(count<=0)
      return(false);

//--- probabilities
   double prob[];
   if(ArrayCopy(prob,probabilities,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate sum and normalize probabilities
   double sum=0;
   for(int i=0; i<size; i++)
      sum+=prob[i];

   if(sum==0.0)
      return(false);

   sum=1.0/sum;

   for(int i=0; i<size; i++)
      prob[i]*=sum;
//--- prepare indices
   int indices[];
   if(ArrayResize(indices,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      indices[i]=i;
//--- sort and prepare cumulative probabilities
   MathQuickSortDescending(prob,indices,0,size-1);
   for(int i=1; i<size; i++)
      prob[i]+=prob[i-1];
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare random samples according to given probabilities
   for(int i=0; i<count; i++)
     {
      double prb=MathRand()/32767.0;
      for(int k=0; k<size; k++)
        {
         if(prb<=prob[k])
           {
            result[i]=array[indices[k]];
            break;
           }
        }
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] according to the given probabilities.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]         : Array with integer values                      |
//| probabilities[] : Array with probabilities                       |
//| count           : Number of elements to choose                   |
//| result[]        : Output array                                   |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const int &array[],double &probabilities[],const int count,int &result[])
  {
   int size=ArraySize(array);
//--- check array size
   if(size==0 || size!=ArraySize(probabilities))
      return(false);

   if(count<=0)
      return(false);

//--- probabilities
   double prob[];
   if(ArrayCopy(prob,probabilities,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate sum and normalize probabilities
   double sum=0;
   for(int i=0; i<size; i++)
      sum+=prob[i];

   if(sum==0.0)
      return(false);

   sum=1.0/sum;

   for(int i=0; i<size; i++)
      prob[i]*=sum;
//--- prepare indices
   int indices[];
   if(ArrayResize(indices,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      indices[i]=i;
//--- sort and prepare cumulative probabilities
   MathQuickSortDescending(prob,indices,0,size-1);
   for(int i=1; i<size; i++)
      prob[i]+=prob[i-1];
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare random samples according to given probabilities
   for(int i=0; i<count; i++)
     {
      double prb=MathRand()/32767.0;
      for(int k=0; k<size; k++)
        {
         if(prb<=prob[k])
           {
            result[i]=array[indices[k]];
            break;
           }
        }
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] according to the given probabilities.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]         : Array with double values                       |
//| probabilities[] : Array with probabilities                       |
//| count           : Number of elements to choose                   |
//| replace         : If true,only unique values are taken to result |
//| result[]        : Output array                                   |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const double &array[],double &probabilities[],const int count,const bool replace,double &result[])
  {
   if(!replace)
      return MathSample(array,probabilities,count,result);

   int size=ArraySize(array);
//--- check array size
   if(size==0 || size!=ArraySize(probabilities))
      return(false);
//--- return false if count>total samples
   if(count<=0 || count>size)
      return(false);

//--- probabilities
   double prob[];
   if(ArrayCopy(prob,probabilities,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate sum and normalize probabilities
   double sum=0;
   for(int i=0; i<size; i++)
      sum+=prob[i];

   if(sum==0.0)
      return(false);

   sum=1.0/sum;

   for(int i=0; i<size; i++)
      prob[i]*=sum;
//--- prepare indices
   int indices[];
   if(ArrayResize(indices,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      indices[i]=i;
//--- sort and prepare cumulative probabilities
   MathQuickSortDescending(prob,indices,0,size-1);
   for(int i=1; i<size; i++)
      prob[i]+=prob[i-1];
//--- prepare array for taken values flag
   bool taken_values[];
   if(ArrayResize(taken_values,size)!=size)
      return(false);
   ArrayFill(taken_values,0,size,false);
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare random samples according to given probabilities
   for(int i=0; i<count; i++)
     {
      double prb=MathRand()/32767.0;
      for(int k=0; k<size; k++)
        {
         int ind=indices[k];
         //--- check also if the value has been already taken
         if(prb<=prob[k] && taken_values[ind]==false)
           {
            result[i]=array[ind];
            taken_values[ind]=true;
            break;
           }
        }
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathSample                                                       |
//+------------------------------------------------------------------+
//| The function takes a sample of the specified size from the       |
//| elements of the array[] according to the given probabilities.    |
//|                                                                  |
//| Arguments:                                                       |
//| array[]         : Array with integer values                      |
//| probabilities[] : Array with probabilities                       |
//| count           : Number of elements to choose                   |
//| replace         : If true,only unique values are taken to result |
//| result[]        : Output array                                   |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathSample(const int &array[],double &probabilities[],const int count,const bool replace,int &result[])
  {
   if(!replace)
      return MathSample(array,probabilities,count,result);

   int size=ArraySize(array);
//--- check array size
   if(size==0 || size!=ArraySize(probabilities))
      return(false);
//--- return false if count>total samples
   if(count<=0 || count>size)
      return(false);

//--- probabilities
   double prob[];
   if(ArrayCopy(prob,probabilities,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- calculate sum and normalize probabilities
   double sum=0;
   for(int i=0; i<size; i++)
      sum+=prob[i];

   if(sum==0.0)
      return(false);

   sum=1.0/sum;

   for(int i=0; i<size; i++)
      prob[i]*=sum;
//--- prepare indices
   int indices[];
   if(ArrayResize(indices,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      indices[i]=i;
//--- sort and prepare cumulative probabilities
   MathQuickSortDescending(prob,indices,0,size-1);
   for(int i=1; i<size; i++)
      prob[i]+=prob[i-1];
//--- prepare array for taken values flag
   bool taken_values[];
   if(ArrayResize(taken_values,size)!=size)
      return(false);
   ArrayFill(taken_values,0,size,false);
//--- prepare target array and calculate values
   if(ArraySize(result)<count)
      if(ArrayResize(result,count)!=count)
         return(false);
//--- prepare random samples according to given probabilities
   for(int i=0; i<count; i++)
     {
      double prb=MathRand()/32767.0;
      for(int k=0; k<size; k++)
        {
         int ind=indices[k];
         //--- check also if the value has been already taken
         if(prb<=prob[k] && taken_values[ind]==false)
           {
            result[i]=array[ind];
            taken_values[ind]=true;
            break;
           }
        }
     }
//---
   return true;
  }
//+------------------------------------------------------------------+
//| MathTukeySummary                                                 |
//+------------------------------------------------------------------+
//| The function calculates Tukey's five number summary.             |
//| It consists of the five most important sample percentiles:       |
//| 1) the sample minimum (smallest observation)                     |
//| 2) the lower quartile or first quartile (0.25 quantile)          |
//| 3) the median (middle value)                                     |
//| 4) the upper quartile or third quartile (0.75 quantile)          |
//| 5) the sample maximum (largest observation)                      |
//|                                                                  |
//| Arguments:                                                       |
//| array[]         : Array with double values                       |
//| removeNAN       : Flag, if true, all NaN values will be removed  |
//| minimum         : Output variable for minimum value              |
//| lower_hinge     : Output variable for 0.25 quantile              |
//| median          : Output variable for median value               |
//| upper_hinge     : Output variable for 0.75 quantile              |
//| maximum         : Output variable for maximum value              |
//|                                                                  |
//| Return value: true if successful, otherwise false.               |
//+------------------------------------------------------------------+
bool MathTukeySummary(const double &array[],const bool removeNAN,double &minimum,double &lower_hinge,double &median,double &upper_hinge,double &maximum)
  {
//--- set default values
   minimum=QNaN;
   lower_hinge=QNaN;
   median=QNaN;
   upper_hinge=QNaN;
   maximum=QNaN;
//--- check array size
   int size=ArraySize(array);
   if(size==0)
      return(false);

   double data[];
   if(ArrayResize(data,size)!=size)
      return(false);

   int actual_values=0;
   if(removeNAN==true)
     {
      //--- remove NAN values (copy non-NAN elements)
      for(int i=0; i<size; i++)
        {
         if(MathIsValidNumber(array[i]))
           {
            data[actual_values]=array[i];
            actual_values++;
           }
        }
      //--- set actual size
      ArrayResize(data,actual_values);
     }
   else
     {
      //--- check NAN and return if NAN found
      for(int i=0; i<size; i++)
        {
         if(!MathIsValidNumber(array[i]))
            return(false);
        }
      //--- otherwise copy data
      if(ArrayCopy(data,array,0,0,WHOLE_ARRAY)!=size)
         return(false);
     }
//--- sort array
   ArraySort(data);
   int n=ArraySize(data);
   if(n==0)
      return(false);
//--- prepare indices
   double n4=MathFloor((n+3)*0.5)*0.5;
   double ind[5];
   ind[0]=0;
   ind[1]=n4-1;
   ind[2]=0.5*(n+1)-1;
   ind[3]=n-n4;
   ind[4]=n-1;
//--- calculate values  
   minimum=0.5*(data[(int)MathFloor(ind[0])]+data[(int)MathCeil(ind[0])]);
   lower_hinge=0.5*(data[(int)MathFloor(ind[1])]+data[(int)MathCeil(ind[1])]);
   median=0.5*(data[(int)MathFloor(ind[2])]+data[(int)MathCeil(ind[2])]);
   upper_hinge=0.5*(data[(int)MathFloor(ind[3])]+data[(int)MathCeil(ind[3])]);
   maximum=0.5*(data[(int)MathFloor(ind[4])]+data[(int)MathCeil(ind[4])]);
//--- successful
   return(true);
  }
//+------------------------------------------------------------------+
//| Computes the minimum and maximum of the values in array[]        |
//+------------------------------------------------------------------+
bool MathRange(const double &array[],double &min,double &max)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- default values, find minimum and maximum values
   min=array[0];
   max=array[0];
   for(int i=1; i<size; i++)
     {
      double value=array[i];
      min=MathMin(min,value);
      max=MathMax(max,value);
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| Computes the minimum value in array[]                            |
//+------------------------------------------------------------------+
double MathMin(const double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(QNaN);
   double min_value=array[0];
   for(int i=1; i<size; i++)
      min_value=MathMin(min_value,array[i]);
//--- return minimum value
   return(min_value);
  }
//+------------------------------------------------------------------+
//| Computes the maximum value in array[]                            |
//+------------------------------------------------------------------+
double MathMax(const double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(QNaN);
   double max_value=array[0];
   for(int i=1; i<size; i++)
      max_value=MathMax(max_value,array[i]);
//--- return maximum value
   return(max_value);
  }
//+------------------------------------------------------------------+
//| Computes the sum of the values in array[]                        |
//+------------------------------------------------------------------+
double MathSum(const double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(QNaN);
//--- calculate sum
   double sum=0.0;
   for(int i=0; i<sum; i++)
      sum+=array[i];
//--- return sum
   return(sum);
  }
//+------------------------------------------------------------------+
//| Computes the product of the values in array[]                    |
//+------------------------------------------------------------------+
double MathProduct(const double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(QNaN);
//--- calculate product
   double product=1.0;
   for(int i=0; i<size; i++)
      product*=array[i];
//--- return product
   return(product);
  }
//+------------------------------------------------------------------+
//| Computes the standard deviation of the values in array[]         |
//+------------------------------------------------------------------+
double MathStandardDeviation(const double &array[])
  {
   int size=ArraySize(array);
   if(size<=1)
      return(QNaN);
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
//--- average mean
   mean=mean/size;
//--- calculate standard deviation   
   double sdev=0;
   for(int i=0; i<size; i++)
      sdev+=MathPow(array[i]-mean,2);
//--- return standard deviation
   return MathSqrt(sdev/(size-1));
  }
//+------------------------------------------------------------------+
//| Computes the average absolute deviation of the values in array[] |
//+------------------------------------------------------------------+
double MathAverageDeviation(const double &array[])
  {
   int size=ArraySize(array);
   if(size<=1)
      return(QNaN);
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
   mean=mean/size;
//--- calculate average deviation
   double adev=0;
   for(int i=0; i<size; i++)
      adev+=MathAbs(array[i]-mean);
   adev=adev/size;
//--- return average deviation
   return(adev);
  }
//+------------------------------------------------------------------+
//| Computes the median value of the values in array[]               |
//+------------------------------------------------------------------+
double MathMedian(double &array[])
  {
   int size=ArraySize(array);
//--- check data range
   if(size==0)
      return(QNaN);
//--- prepare sorted values
   double sorted_values[];
   if(ArrayCopy(sorted_values,array,0,0,WHOLE_ARRAY)!=size)
      return(QNaN);
   ArraySort(sorted_values);
//--- calculate median for odd and even cases
//--- data_count=odd
   if(size%2==1)
      return(sorted_values[size/2]);
   else
//--- data_count=even
      return(0.5*(sorted_values[(size-1)/2]+sorted_values[(size+1)/2]));
  }
//+------------------------------------------------------------------+
//| Computes the mean value of the values in array[]                 |
//+------------------------------------------------------------------+
double MathMean(const double &array[])
  {
   int size=ArraySize(array);
//--- check data range
   if(size<1)
      return(QNaN); // need at least 1 observation
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
   mean=mean/size;
//--- return mean
   return(mean);
  }
//+------------------------------------------------------------------+
//| Computes the variance of the values in array[]                   |
//+------------------------------------------------------------------+
double MathVariance(const double &array[])
  {
   int size=ArraySize(array);
//--- check data range
   if(size<2)
      return(QNaN); // need at least 2 observations
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
   mean=mean/size;
//--- calculate variance
   double variance=0;
   for(int i=0; i<size; i++)
      variance+=MathPow(array[i]-mean,2);
   variance=variance/(size-1);
//--- return variance
   return(variance);
  }
//+------------------------------------------------------------------+
//| Computes the skewness of the values in array[]                   |
//+------------------------------------------------------------------+
double MathSkewness(const double &array[])
  {
   int size=ArraySize(array);
//--- check data range
   if(size<3)
      return(QNaN); // need at least 3 observations
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
   mean=mean/size;
//--- calculate variance and skewness
   double variance=0;
   double skewness=0;
   for(int i=0; i<size; i++)
     {
      double sqr_dev=MathPow(array[i]-mean,2);
      skewness+=sqr_dev*(array[i]-mean);
      variance+=sqr_dev;
     }
   variance=(variance)/(size-1);
   double v3=MathPow(MathSqrt(variance),3);
//---
   if(v3!=0)
     {
      skewness=skewness/(size*v3);
      //--- return skewness
      return(skewness);
     }
   else
      return(QNaN);
  }
//+------------------------------------------------------------------+
//| Computes the kurtosis of the values in array[]                   |
//+------------------------------------------------------------------+
double MathKurtosis(const double &array[])
  {
   int size=ArraySize(array);
//--- check data range
   if(size<4)
      return(QNaN); // need at least 4 observations
//--- calculate mean
   double mean=0.0;
   for(int i=0; i<size; i++)
      mean+=array[i];
   mean=mean/size;
//--- calculate variance and kurtosis
   double variance=0;
   double kurtosis=0;
   for(int i=0; i<size; i++)
     {
      double sqr_dev=MathPow(array[i]-mean,2);
      variance+=sqr_dev;
      kurtosis+=sqr_dev*sqr_dev;
     }
//--- calculate variance
   variance=(variance)/(size-1);
   double v4=MathPow(MathSqrt(variance),4);

   if(v4!=0)
     {
      //--- calculate kurtosis
      kurtosis=kurtosis/(size*v4);
      kurtosis-=3;
      //--- return kurtosis
      return(kurtosis);
     }
   else
      return(QNaN);
  }
//+------------------------------------------------------------------+
//| MathLog1p                                                        |
//+------------------------------------------------------------------+
//| The function calculates log(1+x) for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog1p(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathLog1p(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathLog1p                                                        |
//+------------------------------------------------------------------+
//| The function calculates log(1+x) for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathLog1p(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathLog1p(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathExpm1                                                        |
//+------------------------------------------------------------------+
//| The function calculates exp(x)-1 for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathExpm1(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathExpm1(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathExpm1                                                        |
//+------------------------------------------------------------------+
//| The function calculates exp(x)-1 for the elements from the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathExpm1(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathExpm1(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSinh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic sine for the elements         |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSinh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathSinh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSinh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic sine for the elements         |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSinh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathSinh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCosh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic cosine for the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCosh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathCosh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCosh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic cosine for the elements       |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCosh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathCosh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTanh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic tangent for the elements      |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTanh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathTanh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathTanh                                                         |
//+------------------------------------------------------------------+
//| The function calculates hyperbolic tangent for the elements      |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathTanh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathTanh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArcsinh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic sine for the elements |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArcsinh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArcsinh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArcsinh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic sine for the elements |
//| from the array[].                                                |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArcsinh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArcsinh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArccosh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic cosine for            |
//| the elements from the array[].                                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArccosh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArccosh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArccosh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic cosine for            |
//| the elements from the array[].                                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArccosh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArccosh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArctanh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic tangent for           |
//| the elements from the array[].                                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArctanh(const double &array[],double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathArctanh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathArctanh                                                      |
//+------------------------------------------------------------------+
//| The function calculates inverse hyperbolic tangent for           |
//| the elements from the array[].                                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with values                                  |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathArctanh(double &array[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathArctanh(array[i]);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| Rounds the value to the specified number of significant digits   |
//+------------------------------------------------------------------+
double MathSignif(const double x,const int digits)
  {
   if(x==0)
      return x;

   int dig=digits;
   if(dig>30)
      return x;
   else
   if(dig<1)
      dig=1;

   double sign=1.0;
   double xx=x;
   if(xx<0.0)
     {
      sign=-sign;
      xx=-xx;
     }
//--- calculate log
   double l10 = MathLog10(xx);
   double e10 = (int)(dig-1-MathFloor(l10));
   double value=0,pwr10;
   if(e10>0)
     {
      pwr10=MathPow(10,e10);
      value=MathFloor((xx*pwr10)+0.5)/pwr10;
     }
   else
     {
      pwr10=MathPow(10,-e10);
      value=MathFloor((xx/pwr10)+0.5)*pwr10;
     }
   return(sign*value);
  }
//+------------------------------------------------------------------+
//| MathSignif                                                       |
//+------------------------------------------------------------------+
//| The function calculates the rounded values (to the specified     |
//| number of significant digits) for the elements from the array[]. |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| digits      : Number of significant digits                       |
//| result[]    : Output array with calculated values                |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSignif(const double &array[],int digits,double &result[])
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- prepare target array
   if(ArraySize(result)<size)
      if(ArrayResize(result,size)!=size)
         return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      result[i]=MathSignif(array[i],digits);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathSignif                                                       |
//+------------------------------------------------------------------+
//| The function calculates the rounded values (to the specified     |
//| number of significant digits) for the elements from the array[]. |
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Array with double values                           |
//| digits      : Number of significant digits                       |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathSignif(double &array[],int digits)
  {
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- calculate values
   for(int i=0; i<size; i++)
      array[i]=MathSignif(array[i],digits);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathRank                                                         |
//+------------------------------------------------------------------+
//| The function calculates tied ranks for the elements of the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Input array with double values                     |
//| rank[]      : Output array with ranked values                    |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathRank(const double &array[],double &rank[])
  {
   int size=ArraySize(array);
   if(size<1)
      return(false);
   if(size==1)
     {
      if(ArrayResize(rank,size)!=size)
         return(false);
      rank[0]=1;
      return(true);
     }
//--- prepare arrays
   double values[];
   int indices[];
   if(ArrayCopy(values,array,0,0,WHOLE_ARRAY)!=size)
      return(false);
   if(!MathSequenceByCount(0,size,size,indices))
      return(false);
   if(ArrayResize(rank,size)!=size)
      return(false);
//---
   int i,j,k,t,tmpi;
   double tmp;
//--- sort
   if(size!=1)
     {
      i=2;
      do
        {
         t=i;
         while(t!=1)
           {
            k=t/2;
            if(values[k-1]>=values[t-1])
               t=1;
            else
              {
               //--- swap
               tmp=values[k-1];
               values[k-1]=values[t-1];
               values[t-1]=tmp;
               tmpi=indices[k-1];
               indices[k-1]=indices[t-1];
               indices[t-1]=tmpi;
               t=k;
              }
           }
         i=i+1;
        }
      while(i<=size);
      i=size-1;
      do
        {
         //--- swap
         tmp=values[i];
         values[i]=values[0];
         values[0]=tmp;
         tmpi=indices[i];
         indices[i]=indices[0];
         indices[0]=tmpi;
         t=1;
         while(t!=0)
           {
            k=2*t;
            if(k>i)
               t=0;
            else
              {
               if(k<i)
                  if(values[k]>values[k-1])
                     k++;
               if(values[t-1]>=values[k-1])
                  t=0;
               else
                 {
                  //--- swap
                  tmp=values[k-1];
                  values[k-1]=values[t-1];
                  values[t-1]=tmp;
                  tmpi=indices[k-1];
                  indices[k-1]=indices[t-1];
                  indices[t-1]=tmpi;
                  t=k;
                 }
              }
           }
         i=i-1;
        }
      while(i>=1);
     }
//--- compute tied ranks
   i=0;
   while(i<size)
     {
      j=i+1;
      while(j<size)
        {
         //--- check
         if(values[j]!=values[i])
            break;
         j=j+1;
        }
      for(k=i;k<j;k++)
         values[k]=1+(i+j-1)*0.5;
      i=j;
     }
//--- set output values
   for(i=0; i<size; i++)
      rank[indices[i]]=values[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathRank                                                         |
//+------------------------------------------------------------------+
//| The function calculates tied ranks for the elements of the array.|
//|                                                                  |
//| Arguments:                                                       |
//| array[]     : Input array with integer values                    |
//| rank[]      : Output array with ranked values                    |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathRank(const int &array[],double &rank[])
  {
   int size=ArraySize(array);
   if(size<1)
      return(false);
   if(size==1)
     {
      if(ArrayResize(rank,size)!=size)
         return(false);
      rank[0]=1;
      return(true);
     }
//--- prepare arrays
   double values[];
   int indices[];

   if(ArrayResize(values,size)!=size)
      return(false);
   for(int i=0; i<size; i++)
      values[i]=array[i];

   if(!MathSequenceByCount(0,size,size,indices))
      return(false);
   if(ArrayResize(rank,size)!=size)
      return(false);
//---
   int i,j,k,t,tmpi;
   double tmp;
//--- sort
   if(size!=1)
     {
      i=2;
      do
        {
         t=i;
         while(t!=1)
           {
            k=t/2;
            if(values[k-1]>=values[t-1])
               t=1;
            else
              {
               //--- swap
               tmp=values[k-1];
               values[k-1]=values[t-1];
               values[t-1]=tmp;
               tmpi=indices[k-1];
               indices[k-1]=indices[t-1];
               indices[t-1]=tmpi;
               t=k;
              }
           }
         i=i+1;
        }
      while(i<=size);
      i=size-1;
      do
        {
         //--- swap
         tmp=values[i];
         values[i]=values[0];
         values[0]=tmp;
         tmpi=indices[i];
         indices[i]=indices[0];
         indices[0]=tmpi;
         t=1;
         while(t!=0)
           {
            k=2*t;
            if(k>i)
               t=0;
            else
              {
               if(k<i)
                  if(values[k]>values[k-1])
                     k++;
               if(values[t-1]>=values[k-1])
                  t=0;
               else
                 {
                  //--- swap
                  tmp=values[k-1];
                  values[k-1]=values[t-1];
                  values[t-1]=tmp;
                  tmpi=indices[k-1];
                  indices[k-1]=indices[t-1];
                  indices[t-1]=tmpi;
                  t=k;
                 }
              }
           }
         i=i-1;
        }
      while(i>=1);
     }
//--- compute tied ranks
   i=0;
   while(i<size)
     {
      j=i+1;
      while(j<size)
        {
         //--- check
         if(values[j]!=values[i])
            break;
         j=j+1;
        }
      for(k=i;k<j;k++)
         values[k]=1+(i+j-1)*0.5;
      i=j;
     }
//--- set output values
   for(i=0; i<size; i++)
      rank[indices[i]]=values[i];
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCorrelationPearson                                           |
//+------------------------------------------------------------------+
//| The function calculates Pearson product-moment correlation       |
//| coefficient.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with double values                     |
//| array2[]    : Second array with double values                    |
//| r           : Variable for calculated value of the Pearson r     |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathCorrelationPearson(const double &array1[],const double &array2[],double &r)
  {
   r=QNaN;
   int size=ArraySize(array1);
   if(size<=1 || ArraySize(array2)!=size)
      return(false);
//--- create variables
   double xmean=0;
   double ymean=0;
   double v=1.0/(double)size;
   double x0=array1[0];
   double y0=array2[0];
   double s=0;
   double xv=0;
   double yv=0;
   double t1=0;
   double t2=0;
   bool   samex=true;
   bool   samey=true;
//--- additonally we calculate SameX and SameY - flag variables which are set to true 
//--- when all X[] (or Y[]) contain exactly same value.
//--- if at least one of them is true, we return zero 
//--- (othwerwise we risk to get nonzero correlation because of roundoff).
   for(int i=0; i<size; i++)
     {
      s=array1[i];
      samex=samex && s==x0;
      xmean+=s*v;
      s=array2[i];
      samey=samey && s==y0;
      ymean+=s*v;
     }
//--- check
   if(samex || samey)
      return(false);
//--- calculation
   s=0;
   for(int i=0; i<size; i++)
     {
      t1=array1[i]-xmean;
      t2=array2[i]-ymean;
      xv+=t1*t1;
      yv+=t2*t2;
      s+=t1*t2;
     }
//--- check
   if(xv==0 || yv==0)
      return(false);
//---
   r=s/MathSqrt(xv*yv);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCorrelationPearson                                           |
//+------------------------------------------------------------------+
//| The function calculates Pearson product-moment correlation       |
//| coefficient.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with integer values                    |
//| array2[]    : Second array with integer values                   |
//| r           : Variable for calculated value of the Pearson r     |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathCorrelationPearson(const int &array1[],const int &array2[],double &r)
  {
   r=QNaN;
   int size=ArraySize(array1);
   if(size<=1 || ArraySize(array2)!=size)
      return(false);
//--- create variables
   double xmean=0;
   double ymean=0;
   double v=1.0/(double)size;
   double x0=array1[0];
   double y0=array2[0];
   double s=0;
   double xv=0;
   double yv=0;
   double t1=0;
   double t2=0;
   bool   samex=true;
   bool   samey=true;
//--- additonally we calculate SameX and SameY - flag variables which are set to true 
//--- when all X[] (or Y[]) contain exactly same value.
//--- if at least one of them is true, we return zero 
//--- (othwerwise we risk to get nonzero correlation because of roundoff).
   for(int i=0; i<size; i++)
     {
      s=array1[i];
      samex=samex && s==x0;
      xmean+=s*v;
      s=array2[i];
      samey=samey && s==y0;
      ymean+=s*v;
     }
//--- check
   if(samex || samey)
      return(false);
//--- calculation
   s=0;
   for(int i=0; i<size; i++)
     {
      t1=array1[i]-xmean;
      t2=array2[i]-ymean;
      xv+=t1*t1;
      yv+=t2*t2;
      s+=t1*t2;
     }
//--- check
   if(xv==0 || yv==0)
      return(false);
//---
   r=s/MathSqrt(xv*yv);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCorrelationSpearman                                          |
//+------------------------------------------------------------------+
//| The function calculates Spearman correlation coefficient.        |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with double values                     |
//| array2[]    : Second array with double values                    |
//| r           : Variable for calculated value of the Spearman r    |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathCorrelationSpearman(const double &array1[],const double &array2[],double &r)
  {
   r=QNaN;
   int size=ArraySize(array1);
   if(size<1 || ArraySize(array2)!=size)
      return(false);
//--- calculate ranks
   double rank_x[];
   double rank_y[];
   if(!MathRank(array1,rank_x))
      return(false);
   if(!MathRank(array2,rank_y))
      return(false);
//---
   return MathCorrelationPearson(rank_x,rank_y,r);
  }
//+------------------------------------------------------------------+
//| MathCorrelationSpearman                                          |
//+------------------------------------------------------------------+
//| The function calculates Spearman correlation coefficient.        |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with integer values                    |
//| array2[]    : Second array with integer values                   |
//| r           : Variable for calculated value of the Spearman r    |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
//| Code from ALGLIB numerical analysis and data processing library  |
//+------------------------------------------------------------------+
bool MathCorrelationSpearman(const int &array1[],const int &array2[],double &r)
  {
   r=QNaN;
   int size=ArraySize(array1);
   if(size<1 || ArraySize(array2)!=size)
      return(false);
//--- calculate ranks
   double rank_x[];
   double rank_y[];
   if(!MathRank(array1,rank_x))
      return(false);
   if(!MathRank(array2,rank_y))
      return(false);
//---
   return MathCorrelationPearson(rank_x,rank_y,r);
  }
//+------------------------------------------------------------------+
//| MathCorrelationKendall                                           |
//+------------------------------------------------------------------+
//| The function calculates Kendall's tau correlation coefficient.   |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with double values                     |
//| array2[]    : Second array with double values                    |
//| tau         : Variable for calculated value of the Kendall's tau |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCorrelationKendall(const double &array1[],const double &array2[],double &tau)
  {
   tau=QNaN;
   int size=ArraySize(array1);
   if(size==0 || ArraySize(array2)!=size)
      return(false);
//---
   int cnt1=0,cnt2=0,cnt=0;
//---
   for(int i=0; i<size; i++)
     {
      for(int j=i+1; j<size; j++)
        {
         double delta1=array1[i]-array1[j];
         double delta2=array2[i]-array2[j];
         double delta=delta1*delta2;
         if(delta==0)
           {
            if(delta1!=0)
               cnt1++;
            if(delta2!=0)
               cnt2++;
           }
         else
           {
            cnt1++;
            cnt2++;
            if(delta>0.0)
               cnt++;
            else
               cnt--;
           }
        }
     }
//--- calculate Kendall tau
   double den=cnt1*cnt2;
   if(den==0)
      return(false);
   tau=cnt/MathSqrt(den);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCorrelationKendall                                           |
//+------------------------------------------------------------------+
//| The function calculates Kendall's tau correlation coefficient.   |
//|                                                                  |
//| Arguments:                                                       |
//| array1[]    : First array with integer values                    |
//| array2[]    : Second array with integer values                   |
//| tau         : Variable for calculated value of the Kendall's tau |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCorrelationKendall(const int &array1[],const int &array2[],double &tau)
  {
   tau=QNaN;
   int size=ArraySize(array1);
   if(size==0 || ArraySize(array2)!=size)
      return(false);
//---
   int cnt1=0,cnt2=0,cnt=0;
//---
   for(int i=0; i<size; i++)
     {
      for(int j=i+1; j<size; j++)
        {
         double delta1=array1[i]-array1[j];
         double delta2=array2[i]-array2[j];
         double delta=delta1*delta2;
         if(delta==0)
           {
            if(delta1!=0)
               cnt1++;
            if(delta2!=0)
               cnt2++;
           }
         else
           {
            cnt1++;
            cnt2++;
            if(delta>0.0)
               cnt++;
            else
               cnt--;
           }
        }
     }
//--- calculate Kendall tau
   double den=cnt1*cnt2;
   if(den==0)
      return(false);
   tau=cnt/MathSqrt(den);
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathQuantile                                                     |
//+------------------------------------------------------------------+
//| The function produces sample quantiles corresponding to the      |
//| given probabilities. The smallest observation corresponds        |
//| to a probability of 0.0 and the largest to a probability of 1.0. |
//+------------------------------------------------------------------+
//| The function calculates estimates of underlying distribution     |
//| quantiles based on one or two order statistics from the supplied |
//| elements in x at probabilities in probs.                         |
//|                                                                  |
//| Reference:                                                       |
//| Hyndman, R. J. and Fan, Y. (1996)                                |
//| "Sample quantiles in statistical packages",                      |
//| American Statistician, vol.50, pp. 361–365.                      |
//+------------------------------------------------------------------+
//| All sample quantiles are defined as weighted averages            |
//| of consecutive order statistics.                                 |
//| Sample quantiles of type i are defined by:                       |
//| Q[i](p) = (1 - gamma)*x[j] + gamma*x[j+1]                        |
//|                                                                  |
//| The function MathQuantile calculates quantiles according to      |
//| default type in R (type=7)                                       |
//| p(k) = (k - 1)/(n - 1).                                          |
//| In this case, p(k) = mode[F(x[k])]. This is used by S.           |
//+------------------------------------------------------------------+
bool MathQuantile(const double &array[],const double &probs[],double &quantile[])
  {
   int size=ArraySize(array);
   int size_p=ArraySize(probs);

   if(size==0 || size_p==0)
      return(false);

//--- check probability
   for(int i=0; i<size_p; i++)
     {
      //--- check NaN
      if(!MathIsValidNumber(probs[i]))
         return(false);
      //--- check probability range
      if(probs[i]<0.0 || probs[i]>1.0)
         return(false);
     }
//--- prepare array with values
   double x_values[];
   if(ArrayCopy(x_values,array,0,0,WHOLE_ARRAY)!=size)
      return(false);
//--- prepare output array
   if(ArraySize(quantile)<size_p)
      if(ArrayResize(quantile,size_p)!=size_p)
         return(false);

   double ind[];
   int lo[];
   int hi[];
   ArrayResize(ind,size_p);
   ArrayResize(lo,size_p);
   ArrayResize(hi,size_p);
//--- calculate ranges
   for(int i=0; i<size_p; i++)
     {
      ind[i]=(size-1)*probs[i];
      lo[i]=(int)MathFloor(ind[i]);
      hi[i]=(int)MathCeil(ind[i]);
     }
//--- calculate quantiles
   int indices[];
   MathSequenceByCount(0,size-1,size,indices);
   MathQuickSort(x_values,indices,0,size-1,1);
   for(int i=0; i<size_p; i++)
     {
      quantile[i]=x_values[lo[i]];
      if(ind[i]>lo[i])
        {
         double gamma=(ind[i]-lo[i]);
         quantile[i]=(1.0-gamma)*quantile[i]+gamma*x_values[hi[i]];
        }
     }
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathProbabilityDensityEmpirical                                  |
//+------------------------------------------------------------------+
//| The function calculates the empirical probability density        |
//| function (pdf) for random values from array[].                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]  : Array with random values                              |
//| count    : Data count, total count pairs (x,pdf(x))              |
//| x[]      : Output array for x values                             |
//| pdf[]    : Output array for empirical pdf(x) values              |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathProbabilityDensityEmpirical(const double &array[],const int count,double &x[],double &pdf[])
  {
   if(count<=1)
      return(false);
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NaN values
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]))
         return(false);
     }
//--- prepare output arrays     
   if(ArraySize(x)<count)
      if(ArrayResize(x,count)!=count)
         return(false);
   if(ArraySize(pdf)<count)
      if(ArrayResize(pdf,count)!=count)
         return(false);
//--- search for min,max and range
   double minv=array[0];
   double maxv=array[0];
   for(int i=1; i<size; i++)
     {
      minv=MathMin(minv,array[i]);
      maxv=MathMax(maxv,array[i]);
     }
   double range=maxv-minv;
   if(range==0)
      return(false);
//--- calculate probability density of the empirical distribution
   for(int i=0; i<count; i++)
     {
      x[i]=minv+i*range/(count-1);
      pdf[i]=0;
     }
   for(int i=0; i<size; i++)
     {
      double v=(maxv-array[i])/range;
      int ind=int((v*(count-1)));
      pdf[ind]++;
     }
//--- normalize values
   double sum=0;
   for(int i=0; i<count; i++)
      sum+=pdf[i];
   if(sum==0)
      return(false);
   double coef=1.0/sum;
   for(int i=0; i<count; i++)
      pdf[i]*=coef;
//---
   return(true);
  }
//+------------------------------------------------------------------+
//| MathCumulativeDistributionEmpirical                              |
//+------------------------------------------------------------------+
//| The function calculates the empirical cumulative distribution    |
//| function (cdf) for random values from array[].                   |
//|                                                                  |
//| Arguments:                                                       |
//| array[]  : Array with random values                              |
//| count    : Data count, total count pairs (x,cdf(x))              |
//| x[]      : Output array for x values                             |
//| cdf[]    : Output array for empirical cdf(x) values              |
//|                                                                  |
//| Return value: true if successful, otherwise false                |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionEmpirical(const double &array[],const int count,double &x[],double &cdf[])
  {
   if(count<=1)
      return(false);
   int size=ArraySize(array);
   if(size==0)
      return(false);
//--- check NaN values
   for(int i=0; i<size; i++)
     {
      if(!MathIsValidNumber(array[i]))
         return(false);
     }
//--- prepare output arrays     
   if(ArraySize(x)<count)
      if(ArrayResize(x,count)!=count)
         return(false);
   if(ArraySize(cdf)<count)
      if(ArrayResize(cdf,count)!=count)
         return(false);
//--- search for min,max and range
   double minv=array[0];
   double maxv=array[0];
   for(int i=1; i<size; i++)
     {
      minv=MathMin(minv,array[i]);
      maxv=MathMax(maxv,array[i]);
     }
   double range=maxv-minv;
   if(range==0)
      return(false);
//--- calculate probability density of the empirical distribution
   double pdf[];
   if(ArrayResize(pdf,count)!=count)
      return(false);
   for(int i=0; i<count; i++)
     {
      x[i]=minv+i*range/(count-1);
      pdf[i]=0;
     }
   for(int i=0; i<size; i++)
     {
      double v=(maxv-array[i])/range;
      int ind=int((v*(count-1)));
      pdf[ind]++;
     }
//--- normalize values
   double sum=0;
   for(int i=0; i<count; i++)
      sum+=pdf[i];
   if(sum==0)
      return(false);
   double coef=1.0/sum;
   for(int i=0; i<count; i++)
      pdf[i]*=coef;
//--- calculate cumulative distribution function
   if(!MathCumulativeSum(pdf,cdf))
      return(false);
//---
   return(true);
  }
//+------------------------------------------------------------------+
