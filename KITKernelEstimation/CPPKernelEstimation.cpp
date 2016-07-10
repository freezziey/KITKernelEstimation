//
//  CPPKernelEstimation.cpp
//  Kernel density estimation with a long double precision
//
//  Created by Christian Fries (KIT_1768416) on 12.06.16.
//  Copyright © 2016 Christian Fries. All rights reserved.
//

#include <Rcpp.h>
#include <sstream>
#include <math.h>

using namespace Rcpp;

//Converts numbers to strings for commands
std::string toString(double value){
    std::ostringstream strs;
    strs << value;
    return strs.str();
}


//[[Rcpp::export]]
long double CPP_gamma(long double value){
    return tgammal(value);
}

//[[Rcpp::export]]
long double CPP_loggamma(long double value){
    return lgammal(value);
}

long double gamma_normalization(long double p,long double b){
    
    long double gamma_value=tgammal(p);
    //TODO solve possible accuracy problems for huge numbes
    return powl(b, p)/gamma_value;
    
}

long double gamma_density(double value, double p, double b, long double normalization=0){
    
    if(normalization==0)normalization=gamma_normalization(p, b);
    
    return powl(value, p-1)*powl(M_E, -(b*value))*normalization;
    
}

//[[Rcpp::export]]
NumericVector CPP_gamma_density(NumericVector value, double p, double b){
    
    NumericVector result;
    
    long double normalization=gamma_normalization(p, b);
    
    for(int i=0;i<value.size();i++){
        
        result.insert(i, gamma_density(value[i], p, b,normalization));
        
    }
    
    return result;
    
}


double KDE_rectangular(double xValue,NumericVector X_i,double h=0,double lowerBoundary=-1.7E+308, double upperBoundary=1.7E+308,int exclude=-1){
    
    if(xValue<lowerBoundary||xValue>upperBoundary)return 0;
    
    else{
        
        int sum=0;
        
        bool CV=false;
        
        for(int i=0;i<X_i.size();i++){
            
            if(i!=exclude){
            
            if((X_i[i]-h<lowerBoundary)&&!(xValue>fabsl(X_i[i]-h))){
                
                sum+=2;
                
            }
            else if(X_i[i]+h>upperBoundary&&!(xValue<fabsl(upperBoundary-(X_i[i]+h-upperBoundary)))) {
                
                sum+=2;
                
            }
            else{
                
                if(!(xValue<X_i[i]-h)&&!(xValue>X_i[i]+h))sum++;
                
            }
        
            }
            
            else CV=true;
            
        }
        
        double result;
        
        (CV) ? result = sum/h/2.0/(X_i.size()-1) : result = sum/h/2.0/X_i.size();
        
        return result;
        
    }
}


//[[Rcpp::export]]
NumericVector CPP_KDE_rectanuglar(NumericVector xValue, NumericVector X_i, double h=0,
                                  double lowerBoundary=-1.7E+308, double upperBoundary=1.7E+308){
    
    NumericVector result;
    
    if(h==0)stop("Please select a bandwith h!");    //TODO automatic bandwith calculation if none is selected
    
    for (int x=0; x<xValue.size(); x++) {
        
        checkUserInterrupt();
        
        result.insert(x, KDE_rectangular(xValue[x], X_i,h,lowerBoundary,upperBoundary));
        
    }
    
    return result;
    
}


double KDE_gamma_MOD(double xValue, NumericVector X_i, double b, int exclude=-1){
    
    long double sum=0.0;
    
    bool CV=false;
    
    for(int i=0;i<X_i.size();i++){
        
        if(X_i[i]<0)stop("Negative value occured in line ["+toString(i)+"]!");
        
        if(i!=exclude){
        
            long double p=X_i[i]*b+1.0;
        
        sum+=gamma_density(xValue, p, b);
            
        }
        
        else CV=true;
        
    }
    
    double result;
    
    (CV) ? result = sum/(X_i.size()-1) : result = sum/X_i.size();
    
    return result;
    
}



//[[Rcpp::export]]
NumericVector CPP_KDE_gamma_MOD(NumericVector xValue, NumericVector X_i, double b=1){
    
    NumericVector result;
    
    for (int x=0; x<xValue.size(); x++) {
        
        checkUserInterrupt();
        
        result.insert(x, KDE_gamma_MOD(xValue[x],X_i,b));
        
    }
    
    return result;
    
}


double KDE_gamma_EXP(double xValue, NumericVector X_i, double b, int exclude=-1){
    
    long double sum=0.0;
    
    bool CV=false;
    
    for(int i=0;i<X_i.size();i++){
        
        if(X_i[i]<0)stop("Negative value occured in line ["+toString(i)+"]!");
        
        if(i!=exclude){
        
        long double p;
        
        (X_i[i]<2/b)? p=powl(X_i[i]*b,2.0)/4.0+1 : p=X_i[i]*b; //Chen(2000) approach
        
        sum+=gamma_density(xValue, p, b);
            
        }
        
        else CV=true;
        
    }
    
    double result;
    
    (CV) ? result = sum/(X_i.size()-1) : result = sum/X_i.size();
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_KDE_gamma_EXP(NumericVector xValue, NumericVector X_i, double b=1){
    
    NumericVector result;
    
    for (int x=0; x<xValue.size(); x++) {
        
        checkUserInterrupt();
        
        result.insert(x, KDE_gamma_EXP(xValue[x],X_i,b));
        
    }
    
    return result;
    
}


double KDE_gammaChen_MOD(double xValue, NumericVector X_i, double b, int exclude=-1){
    
    long double sum=0.0;
    
    long double p=xValue*b+1.0;
    
    long double normalization=gamma_normalization(p, b);
    
    bool CV=false;
    
    for(int i=0;i<X_i.size();i++){
        
        if(X_i[i]<0)Rcpp::stop("Negative value occured in line ["+toString(i)+"]!");
        
        if(i!=exclude){
        
        sum+=gamma_density(X_i[i], p, b,normalization);
            
        }
        
        else CV=true;
        
    }
    
    double result;
    
    (CV) ? result = sum/(X_i.size()-1) : result = sum/X_i.size();
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_KDE_gammaChen_MOD(NumericVector xValue, NumericVector X_i, double b=1){
    
    NumericVector result;
    
    for (int x=0; x<xValue.size(); x++) {
        
        checkUserInterrupt();
        
        result.insert(x, KDE_gammaChen_MOD(xValue[x],X_i,b));
        
    }
    
    return result;
    
}

double KDE_gammaChen_EXP(double xValue, NumericVector X_i, double b, int exclude=-1){
    
    long double sum=0.0;
    
    long double p;
    
    (xValue<2/b)? p=powl(xValue*b,2.0)/4.0+1 : p=xValue*b; //Chen(2000) approach
    
    long double normalization=gamma_normalization(p, b);
    
    bool CV=false;
    
    for(int i=0;i<X_i.size();i++){
        
        if(X_i[i]<0)Rcpp::stop("Negative value occured in line ["+toString(i)+"]!");
        
        if(i!=exclude){
        
        sum+=gamma_density(X_i[i], p, b,normalization);
            
        }
        
        else CV=true;
        
    }
    
    double result;
    
    (CV) ? result = sum/(X_i.size()-1) : result = sum/X_i.size();
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_KDE_gammaChen_EXP(NumericVector xValue, NumericVector X_i, double b=1){
    
    NumericVector result;
    
    for (int x=0; x<xValue.size(); x++) {
        
        checkUserInterrupt();
        
        result.insert(x, KDE_gammaChen_EXP(xValue[x],X_i,b));
        
    }
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_CrossValidation_rectangular(NumericVector bandwidth, NumericVector xValues, NumericVector X_i,double lowerBoundary=-1.7E+308,double upperBoundary=1.7E+308){
    
    NumericVector result;
    
    for (int h=0; h<bandwidth.size(); h++) {
        
        checkUserInterrupt();
        
        if(!(bandwidth[h]>0))stop("Select bandwidth > 0!");
        
        //calculate the integrated_squared_estimation
        long double integrated__squared_estimation=0.0;
            
        NumericVector KDE=CPP_KDE_rectanuglar(xValues, X_i,bandwidth[h],lowerBoundary,upperBoundary);
        
         for (int x=0; x<xValues.size(); x++) {
             
             long double squared_est=powl(KDE[x], 2.0);
             
             integrated__squared_estimation+=squared_est;
             
         }
            
        double area=xValues[xValues.size()-1]-xValues[0];
        
        integrated__squared_estimation=integrated__squared_estimation/xValues.size()*area;
        
        //calculate the cross validation
        long double CV=0.0;
        
        for (int i=0; i<X_i.size(); i++) {
            
            CV+=KDE_rectangular(X_i[i], X_i,bandwidth[h],lowerBoundary,upperBoundary,i);
            
        }
        
        CV=CV/X_i.size()*2.0;
        
        result.insert(h, integrated__squared_estimation-CV);
        
    }
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_CrossValidation_gamma_MOD(NumericVector bandwidth, NumericVector xValues, NumericVector X_i){
    
    NumericVector result;
    
    for (int h=0; h<bandwidth.size(); h++) {
        
        checkUserInterrupt();
        
        if(!(bandwidth[h]>0))stop("Select bandwidth > 0!");
        
        //calculate the integrated_squared_estimation
        long double integrated__squared_estimation=0.0;
        
        NumericVector KDE=CPP_KDE_gamma_MOD(xValues, X_i,bandwidth[h]);
        
        for (int x=0; x<xValues.size(); x++) {
            
            long double squared_est=powl(KDE[x], 2.0);
            
            integrated__squared_estimation+=squared_est;
            
        }
        
        double area=xValues[xValues.size()-1]-xValues[0];
        
        integrated__squared_estimation=integrated__squared_estimation/xValues.size()*area;
        
        //calculate the cross validation
        long double CV=0.0;
        
        for (int i=0; i<X_i.size(); i++) {
            
            CV+=KDE_gamma_MOD(X_i[i], X_i, bandwidth[h],i);
            
        }
        
        CV=CV/X_i.size()*2.0;
        
        
        result.insert(h, integrated__squared_estimation-CV);
        
        
    }
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_CrossValidation_gamma_EXP(NumericVector bandwidth, NumericVector xValues, NumericVector X_i){
    
    NumericVector result;
    
    for (int h=0; h<bandwidth.size(); h++) {
        
        checkUserInterrupt();
        
        if(!(bandwidth[h]>0))stop("Select bandwidth > 0!");
        
        //calculate the integrated_squared_estimation
        long double integrated__squared_estimation=0.0;
        
        NumericVector KDE=CPP_KDE_gamma_EXP(xValues, X_i,bandwidth[h]);
        
        for (int x=0; x<xValues.size(); x++) {
            
            long double squared_est=powl(KDE[x], 2.0);
            
            integrated__squared_estimation+=squared_est;
            
        }
        
        double area=xValues[xValues.size()-1]-xValues[0];
        
        integrated__squared_estimation=integrated__squared_estimation/xValues.size()*area;
        
        //calculate the cross validation
        long double CV=0.0;
        
        for (int i=0; i<X_i.size(); i++) {
            
            CV+=KDE_gamma_EXP(X_i[i], X_i, bandwidth[h],i);
            
        }
        
        CV=CV/X_i.size()*2.0;
        
        
        result.insert(h, integrated__squared_estimation-CV);
        
        
    }
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_CrossValidation_gammaChen_MOD(NumericVector bandwidth, NumericVector xValues, NumericVector X_i){
    
    NumericVector result;
    
    for (int h=0; h<bandwidth.size(); h++) {
        
        checkUserInterrupt();
        
        if(!(bandwidth[h]>0))stop("Select bandwidth > 0!");
        
        //calculate the integrated_squared_estimation
        long double integrated__squared_estimation=0.0;
        
        NumericVector KDE=CPP_KDE_gammaChen_MOD(xValues, X_i,bandwidth[h]);
        
        for (int x=0; x<xValues.size(); x++) {
            
            long double squared_est=powl(KDE[x], 2.0);
            
            integrated__squared_estimation+=squared_est;
            
        }
        
        double area=xValues[xValues.size()-1]-xValues[0];
        
        integrated__squared_estimation=integrated__squared_estimation/xValues.size()*area;
        
        //calculate the cross validation
        long double CV=0.0;
        
        for (int i=0; i<X_i.size(); i++) {
            
            CV+=KDE_gammaChen_MOD(X_i[i], X_i, bandwidth[h],i);
            
        }
        
        CV=CV/X_i.size()*2.0;
        
        
        result.insert(h, integrated__squared_estimation-CV);
        
        
    }
    
    return result;
    
}


//[[Rcpp::export]]
NumericVector CPP_CrossValidation_gammaChen_EXP(NumericVector bandwidth, NumericVector xValues, NumericVector X_i){
    
    NumericVector result;
    
    for (int h=0; h<bandwidth.size(); h++) {
        
        checkUserInterrupt();
        
        if(!(bandwidth[h]>0))stop("Select bandwidth > 0!");
        
        //calculate the integrated_squared_estimation
        long double integrated__squared_estimation=0.0;
        
        NumericVector KDE=CPP_KDE_gammaChen_EXP(xValues, X_i,bandwidth[h]);
        
        for (int x=0; x<xValues.size(); x++) {
            
            long double squared_est=powl(KDE[x], 2.0);
            
            integrated__squared_estimation+=squared_est;
            
        }
        
        double area=xValues[xValues.size()-1]-xValues[0];
        
        integrated__squared_estimation=integrated__squared_estimation/xValues.size()*area;
        
        //calculate the cross validation
        long double CV=0.0;
        
        for (int i=0; i<X_i.size(); i++) {
            
            CV+=KDE_gammaChen_EXP(X_i[i], X_i, bandwidth[h],i);
            
        }
        
        CV=CV/X_i.size()*2.0;
        
        
        result.insert(h, integrated__squared_estimation-CV);
        
    }
    
    return result;
    
}




