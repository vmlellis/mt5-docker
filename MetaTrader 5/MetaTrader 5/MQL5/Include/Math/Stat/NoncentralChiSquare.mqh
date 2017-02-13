//+------------------------------------------------------------------+
//|                                          NoncentralChiSquare.mqh |
//|                        Copyright 2016, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2016, MetaQuotes Software Corp."
#property link      "https://www.mql5.com"

#include "Math.mqh"
#include "Normal.mqh"
#include "Poisson.mqh"
#include "ChiSquare.mqh"
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution density function (PDF)        |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| Noncentral Chi-Square distribution with parameters nu and sigma. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralChiSquare(const double x,const double nu,const double sigma,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(true,log_mode);
//--- factors
   double lambda=sigma*0.5;
   double coef_lambda_x=MathExp(-lambda-x*0.5);
   double half_nu=nu*0.5;
   double pwr_lambda=1.0;
   double pwr_x=MathExp((half_nu-1.0)*MathLog(x));
   double pwr_two=MathExp(-half_nu*MathLog(2));
   double fact_mult=1.0;
   double pwr_gamma=1.0/MathGamma(half_nu);
//--- calculate density using direct summation
   const int max_terms=200;
   int j=0;
   double pdf=0.0;
   while(j<max_terms)
     {
      if(j>0)
        {
         pwr_lambda*=lambda;
         pwr_x*=x;
         pwr_two*=0.5;
         double inv_j=1.0/j;
         fact_mult*=inv_j;
         pwr_gamma*=inv_j;
        }
      double dp=coef_lambda_x*pwr_lambda*pwr_x*pwr_two*pwr_gamma*fact_mult;
      pdf=pdf+dp;
      //--- check stop
      if(dp/(pdf+10E-10)<10E-16)
         break;
      j++;
     }
//--- check convergence
   if(j<max_terms)
      return TailLogValue(pdf,true,log_mode);
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution density function (PDF)        |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| Noncentral Chi-Square distribution with parameters nu and sigma. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralChiSquare(double x,const double nu,const double sigma,int &error_code)
  {
   return MathProbabilityDensityNoncentralChiSquare(x,nu,sigma,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution density function (PDF)        |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Chi Square distribution with parameter nu for values in x.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityNoncentralChiSquare(const double &x[],const double nu,const double sigma,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

//--- common factors
   double lambda=sigma*0.5;
   double half_nu=nu*0.5;
   double pwr_two0=MathExp(-half_nu*MathLog(2));
   double pwr_gamma0=1.0/MathGamma(half_nu);
   double coef_lambda=MathExp(-lambda);
   double half_nu_m1=half_nu-1.0;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         double coef_lambda_x=coef_lambda*MathExp(-x_arg*0.5);
         double pwr_lambda=1;
         double pwr_x=MathExp(half_nu_m1*MathLog(x_arg));
         double fact_mult=1.0;
         double pwr_two=pwr_two0;
         double pwr_gamma=pwr_gamma0;
         //--- calculate density using direct summation
         int max_terms=200;
         int j=0;
         double pdf=0.0;
         while(j<max_terms)
           {
            if(j>0)
              {
               pwr_lambda*=lambda;
               pwr_x*=x_arg;
               pwr_two*=0.5;
               double inv_j=1.0/j;
               fact_mult*=inv_j;
               pwr_gamma*=inv_j;
              }
            double dp=coef_lambda_x*pwr_lambda*pwr_x*pwr_two*pwr_gamma*fact_mult;
            pdf=pdf+dp;
            //--- check stop
            if(dp/(pdf+10E-10)<10E-16)
               break;
            j++;
           }
         //--- check convergence
         if(j<max_terms)
            result[i]=TailLogValue(pdf,true,log_mode);
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution density function (PDF)        |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Chi-Square distribution with parameter nu for values in x[]. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityNoncentralChiSquare(const double &x[],const double nu,const double sigma,double &result[])
  {
   return MathProbabilityDensityNoncentralChiSquare(x,nu,sigma,false,result);
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square cumulative distribution function (CDF)     |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from a Noncentral Chi-Square distribution with parameters        |
//| nu and sigma is less than or equal to x.                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of Noncentral Chi-Square cumulative distribution       |
//| function with parameters nu and sigma, evaluated at x.           |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralChiSquare(const double x,const double nu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(true,log_mode);

//--- prepare parameters
   double cdf=0.0;
   int max_terms=100;

   double lambda=sigma*0.5;
   double coef_lambda=MathExp(-lambda);
   double pwr_lambda=1.0;
   double fact_mult=1.0;
   double half_x=x*0.5;
   double half_nu=nu*0.5;
//--- direct summation
   int j=0;
   while(j<max_terms)
     {
      if(j>0)
        {
         pwr_lambda*=lambda;
         fact_mult/=j;
        }
      double coef1=coef_lambda*pwr_lambda*fact_mult;
      double coef2=MathMin(MathGammaIncomplete(half_x,half_nu+j),1.0);
      double dp=coef1*coef2;
      cdf=cdf+dp;
      if((dp/(cdf+10E-10))<10E-16)
         break;
      j++;
     }
//---
   if(j<max_terms)
     {
      //--- take into account round-off errors for probability
      cdf=MathMin(cdf,1.0);
      return TailLogValue(cdf,tail,log_mode);
     }
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square cumulative distribution function (CDF)     |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from a Noncentral Chi-Square distribution with parameters        |
//| nu and sigma is less than or equal to x.                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of Noncentral Chi-Square cumulative distribution       |
//| function with parameters nu and sigma, evaluated at x.           |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralChiSquare(const double x,const double nu,const double sigma,int &error_code)
  {
   return MathCumulativeDistributionNoncentralChiSquare(x,nu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square cumulative distribution function (CDF)     |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Chi-Square distribution with parameters nu and    |
//| sigma for values in x.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralChiSquare(const double &x[],const double nu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

//--- common factors
   double lambda=sigma*0.5;
   double coef_lambda=MathExp(-lambda);
   double half_nu=nu*0.5;
   const int max_terms=100;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0.0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         double pwr_lambda=1.0;
         double fact_mult=1.0;
         double half_x=x_arg*0.5;
         double cdf=0.0;
         int j=0;
         //--- direct summation
         while(j<max_terms)
           {
            if(j>0)
              {
               pwr_lambda*=lambda;
               fact_mult/=j;
              }
            double coef1=coef_lambda*pwr_lambda*fact_mult;
            double coef2=MathMin(MathGammaIncomplete(half_x,half_nu+j),1.0);
            double dp=coef1*coef2;
            cdf=cdf+dp;
            if((dp/(cdf+10E-10))<10E-16)
               break;
            j++;
           }
         //---
         if(j<max_terms)
           {
            //--- take into account round-off errors for probability
            cdf=MathMin(cdf,1.0);
            result[i]=TailLogValue(cdf,tail,log_mode);
           }
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square cumulative distribution function (CDF)     |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Chi-Square distribution with parameters nu and    |
//| sigma for values in x.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralChiSquare(const double &x[],const double nu,const double sigma,double &result[])
  {
   return MathCumulativeDistributionNoncentralChiSquare(x,nu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution quantile function(inverse CDF)|
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral Chi-Square distribution with parameters   |
//| nu and sigma for the desired probability.                        |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| sigma       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function of     |
//| Noncentral Chi-Square distribution with parameters nu and sigma. |
//+------------------------------------------------------------------+
double MathQuantileNoncentralChiSquare(const double probability,const double nu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
   if(log_mode==true)
     {
      if(probability==QNEGINF)
         return 0.0;
     }
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

//--- calculate real probability
   double prob=TailLogProbability(probability,tail,log_mode);
//--- check probability range
   if(prob<0.0 || prob>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(prob==0.0)
      return 0.0;

   if(prob==1.0)
      return QPOSINF;

   error_code=ERR_OK;

//--- common factors for pdf and cdf calculation
   double lambda=sigma*0.5;
   double half_nu=nu*0.5;
   double pwr_two0=MathExp(-half_nu*MathLog(2));
   double pwr_gamma0=1.0/MathGamma(half_nu);
   double coef_lambda=MathExp(-lambda);
   double half_nu_m1=half_nu-1.0;
   const int max_terms=200;
//--- prepare values for initial x estimation
   double x=0.5;
   double h=1.0;
   double h_min=10E-10;
//--- Newton iterations
   const int max_iterations=50;
   int iterations=0;
   while(iterations<max_iterations)
     {
      //--- check convergence
      if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
         break;

      //--- calculate pdf
      double half_x=x*0.5;
      double coef_lambda_x=coef_lambda*MathExp(-half_x);
      double pwr_lambda=1.0;
      double pwr_x=MathExp(half_nu_m1*MathLog(x));
      double fact_mult=1.0;
      double pwr_two=pwr_two0;
      double pwr_gamma=pwr_gamma0;
      //--- calculate density using direct summation
      int j=0;
      double pdf=0.0;
      while(j<max_terms)
        {
         if(j>0)
           {
            pwr_lambda*=lambda;
            pwr_x*=x;
            pwr_two*=0.5;
            double inv_j=1.0/j;
            fact_mult*=inv_j;
            pwr_gamma*=inv_j;
           }
         double dp=coef_lambda_x*pwr_lambda*pwr_x*pwr_two*pwr_gamma*fact_mult;
         pdf=pdf+dp;
         //--- check stop
         if(dp/(pdf+10E-10)<10E-16)
            break;
         j++;
        }
      //--- check convergence
      if(j>max_terms)
        {
         error_code=ERR_NON_CONVERGENCE;
         return QNaN;
        }

      //--- calculate cdf
      pwr_lambda=1.0;
      fact_mult=1.0;
      double cdf=0.0;
      j=0;
      //--- direct summation
      while(j<max_terms)
        {
         if(j>0)
           {
            pwr_lambda*=lambda;
            fact_mult/=j;
           }
         double coef1=coef_lambda*pwr_lambda*fact_mult;
         double coef2=MathMin(MathGammaIncomplete(half_x,half_nu+j),1.0);
         double dp=coef1*coef2;
         cdf=cdf+dp;
         if((dp/(cdf+10E-10))<10E-16)
            break;
         j++;
        }
      //---
      if(j>max_terms)
        {
         error_code=ERR_NON_CONVERGENCE;
         return QNaN;
        }

      //--- calculate ratio
      h=(cdf-prob)/pdf;

      double x_new=x-h;
      if(x_new<0.0)
         x_new=x*0.1;
      else
      if(x_new>1.0)
         x_new=1.0-(1-x)*0.1;
      x=x_new;

      iterations++;
     }
//--- check convergence
   if(iterations<max_iterations)
      return x;
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution quantile function(inverse CDF)|
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Noncentral Chi-Square distribution               |
//| with parameters mu and sigma for the desired probability.        |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| sigma       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function of     |
//| Noncentral Chi-Square distribution with parameters mu and sigma. |
//+------------------------------------------------------------------+
double MathQuantileNoncentralChiSquare(const double probability,const double nu,const double sigma,int &error_code)
  {
   return MathQuantileNoncentralChiSquare(probability,nu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution quantile function(inverse CDF)|
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral Chi-Square distribution with          |
//| parameters nu and sigma for values from the probability[] array. |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| sigma       : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralChiSquare(const double &probability[],const double nu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

//--- common factors for pdf and cdf calculation
   double lambda=sigma*0.5;
   double half_nu=nu*0.5;
   double pwr_two0=MathExp(-half_nu*MathLog(2));
   double pwr_gamma0=1.0/MathGamma(half_nu);
   double coef_lambda=MathExp(-lambda);
   double half_nu_m1=half_nu-1.0;
   const int max_terms=200;
   const int max_iterations=50;
   double h_min=10E-10;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         //--- check probability range
         if(prob<0.0 || prob>1.0)
            return false;

         //--- prepare values for initial x estimation
         int    err_code=0;
         double x=0.5;
         double h=1.0;
         //--- Newton iterations
         int iterations=0;
         while(iterations<max_iterations)
           {
            //--- check convergence
            if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
               break;

            //--- calculate pdf
            double half_x=x*0.5;
            double coef_lambda_x=coef_lambda*MathExp(-half_x);
            double pwr_lambda=1;
            double pwr_x=MathExp(half_nu_m1*MathLog(x));
            double fact_mult=1.0;
            double pwr_two=pwr_two0;
            double pwr_gamma=pwr_gamma0;
            //--- calculate density using direct summation
            int j=0;
            double pdf=0.0;
            while(j<max_terms)
              {
               if(j>0)
                 {
                  pwr_lambda*=lambda;
                  pwr_x*=x;
                  pwr_two*=0.5;
                  double inv_j=1.0/j;
                  fact_mult*=inv_j;
                  pwr_gamma*=inv_j;
                 }
               double dp=coef_lambda_x*pwr_lambda*pwr_x*pwr_two*pwr_gamma*fact_mult;
               pdf=pdf+dp;
               //--- check stop
               if(dp/(pdf+10E-10)<10E-16)
                  break;
               j++;
              }
            //--- check convergence
            if(j>max_terms)
               return false;

            //--- calculate cdf
            pwr_lambda=1.0;
            fact_mult=1.0;
            double cdf=0.0;
            j=0;
            //--- direct summation
            while(j<max_terms)
              {
               if(j>0)
                 {
                  pwr_lambda*=lambda;
                  fact_mult/=j;
                 }
               double coef1=coef_lambda*pwr_lambda*fact_mult;
               double coef2=MathMin(MathGammaIncomplete(half_x,half_nu+j),1.0);
               double dp=coef1*coef2;
               cdf=cdf+dp;
               if((dp/(cdf+10E-10))<10E-16)
                  break;
               j++;
              }
            //---
            if(j>max_terms)
               return false;

            //--- calculate ratio
            h=(cdf-prob)/pdf;

            double x_new=x-h;
            if(x_new<0.0)
               x_new=x*0.1;
            else
            if(x_new>1.0)
               x_new=1.0-(1-x)*0.1;
            x=x_new;

            iterations++;
           }
         //--- check convergence
         if(iterations<max_iterations)
            result[i]=x;
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution quantile function(inverse CDF)|
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral Chi-Square distribution with          |
//| parameters nu and sigma for values from the probability[] array. |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| sigma       : Noncentrality parameter                            |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralChiSquare(const double &probability[],const double nu,const double sigma,double &result[])
  {
   return MathQuantileNoncentralChiSquare(probability,nu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral Chi-Square distribution       |
//+------------------------------------------------------------------+
//| Compute the random variable from the Noncentral Chi-Square       |
//| distribution with parameters nu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| nu          : Degrees of freedom                                 |
//| sigma       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Noncentral Chi-Square distribution.        |
//+------------------------------------------------------------------+
//| Author: Robert Kern                                              |
//+------------------------------------------------------------------+
double MathRandomNoncentralChiSquare(const double nu,const double sigma,int &error_code)
  {
//--- return ChiSquare if sigma==0
   if(sigma==0.0)
     {
      return MathRandomChiSquare(nu,error_code);
     }
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check sigma
   if(sigma<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(nu>1)
     {
      const double chi2=2.0*MathRandomGamma((nu-1)*0.5,1.0,error_code);
      double f,x1,x2,r2;
      do
        {
         x1=2.0*MathRandomNonZero()-1.0;
         x2=2.0*MathRandomNonZero()-1.0;
         r2=x1*x1+x2*x2;
        }
      while(r2>=1.0 || r2==0.0);
      //--- Box-Muller transform
      f=MathSqrt(-2.0*MathLog(r2)/r2);
      const double N2=MathPow(f*x2+MathSqrt(sigma),2);
      return chi2+N2;
     }
   else
     {
      double x=MathRandomPoisson(sigma*0.5,error_code);
      return 2.0*MathRandomGamma((nu+2*x)*0.5,1.0,error_code);
     }
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral Chi-Square distribution       |
//+------------------------------------------------------------------+
//| Generates random variables from the Noncentral Chi-Square        |
//| distribution with parameters nu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNoncentralChiSquare(const double nu,const double sigma,const int data_count,double &result[])
  {
//--- return ChiSquare if sigma==0
   if(sigma==0.0)
      return MathRandomChiSquare(nu,data_count,result);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
      return false;
//--- check nu
   if(nu!=MathRound(nu) || nu<=0)
      return false;
//--- check sigma
   if(sigma<0.0)
      return false;

   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      if(nu>1)
        {
         const double chi2=2.0*MathRandomGamma((nu-1)*0.5,1.0,error_code);
         double f,x1,x2,r2;
         do
           {
            x1=2.0*MathRandomNonZero()-1.0;
            x2=2.0*MathRandomNonZero()-1.0;
            r2=x1*x1+x2*x2;
           }
         while(r2>=1.0 || r2==0.0);
         //--- Box-Muller transform
         f=MathSqrt(-2.0*MathLog(r2)/r2);
         const double N2=MathPow(f*x2+MathSqrt(sigma),2);
         result[i]=chi2+N2;
        }
      else
        {
         double x=MathRandomPoisson(sigma*0.5,error_code);
         result[i]=2.0*MathRandomGamma((nu+2*x)*0.5,1.0,error_code);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral Chi-Square distribution moments                       |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Noncental Chi-Square  |
//| distribution with parameters nu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| sigma      : Noncentrality parameter                             |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsNoncentralChiSquare(const double nu,const double sigma,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =nu+sigma;
   variance=2*nu+4*sigma;
   skewness=2*M_SQRT2*(nu+3*sigma)*MathPow(nu+2*sigma,-1.5);
   kurtosis=12*(nu+4*sigma)*MathPow(nu+2*sigma,-2);
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
