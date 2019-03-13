/*
 *Copyright 2019 Anders Reenberg Andersen, PhD
 *
 *Licensed under the Apache License, Version 2.0 (the "License");
 *you may not use this file except in compliance with the License.
 *You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 *distributed under the License is distributed on an "AS IS" BASIS,
 *WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *See the License for the specific language governing permissions and
 *limitations under the License.
 */

package queueing;
/*
 * This class evaluates the state probability distribution of  
 * a time-dependent M/M/C/K queueing system with homogeneous parameter input.
 * Author contact details: andersra@live.dk, areenbergandersen@gmail.com
 *
 *
 * @author Anders Reenberg Andersen, PhD
 */
public class evaluate {
    
    
  int Ns; //total number of states
  double gamma; //uniformization rate
  double lambda; //the constant arrival rate
  double mu; //the constant service rate
  int servers; //the number of servers
  double[][] Q;
  double[] pi0; //initial state probability distribution
  double[] pi; //state probability distribution after time t
  int currentlyUsed; //current number of non-idle servers
  
  
  
  public evaluate(int ch,double l, double m, int c, int cap){ //class constructor with non-idle servers as input 
      
      currentlyUsed = ch;
      lambda = l;
      mu = m;
      servers = c;
      
      Ns = cap + 1;
      
      transitionratematrix();
      stochasticmatrix();
      initialDistribution();
  }
  
  
  public evaluate(double[] initDist,double l, double m, int c, int cap){ //class constructor with distribution input
      
      lambda = l;
      mu = m;
      servers = c;
      
      Ns = cap + 1;
      
      transitionratematrix();
      stochasticmatrix();
      
      pi0 = new double[Ns];
      int i;
      
      if (Ns>initDist.length){
        for(i=0; i<initDist.length; i++){
            pi0[i] = initDist[i];
        }
        for(i=initDist.length; i<Ns; i++){
            pi0[i] = 0;
        }
      }else if(Ns<initDist.length){
          for(i=0; i<=(Ns-2); i++){
              pi0[i] = initDist[i];
          }
          double sm = 0;
          for(i=(Ns-1); i<initDist.length; i++){
              sm += initDist[i];
          }
          pi0[(Ns-1)] = sm;
      }else{
        for(i=0; i<Ns; i++){
            pi0[i] = initDist[i];
        }  
      }
      
      //System.out.println("initDist");
      //for(i=0; i<initDist.length; i++){
      //    System.out.print(initDist[i] + " ");
      //}
      //System.out.println();
      
      
      //System.out.println("pi0");
      //for(i=0; i<pi0.length; i++){
      //    System.out.print(pi0[i] + " ");
      //}
      //System.out.println();
  }
  
  
  private void initialDistribution(){
      
      pi0 = new double[Ns];
      int i;
      for(i=0; i<Ns; i++){
          
          if (i == (Ns-1) && currentlyUsed >= Ns){
              pi0[i] = 1;
          }else if (i == currentlyUsed && currentlyUsed < Ns){
              pi0[i] = 1;
          }else{
              pi0[i] = 0;
          }
          
      }
      
  }
  
  
    
  private void transitionratematrix(){
        //generates the transition rate matrix for a M/M/C/K queueing system in Harwell-Boeing format
        
      
        int nz = Ns+(Ns-1)*2;
        gamma = 0;

        int i; //generate state space
        int[] S = new int[Ns];
        for(i=0; i<Ns; i++){
            S[i] = i;
        }
        
                
        Q = new double[3][];
        Q[0] = new double[nz]; //aa
        Q[1] = new double[nz]; //ja
        Q[2] = new double[Ns+1]; //ia
        
        //----------------------------------------
        
        
        int v; //state to be evaluated
        int rn = 1; Q[2][0] = rn; double dv; double val;
        gamma = 0;
        for(i=0; i<Ns; i++){
            
            //get state
            v = S[i];
            dv = 0;
            
            //evaluate each jump-type 
            if(v<(Ns-1)){ //forward
                val = lambda;
                
                dv += val;
                Q[0][rn-1] = val;
                Q[1][rn-1] = (i+1)+1;
                rn += 1;        
            }
            if(v>0){ //backward
                if (v<=servers){
                    val = mu*v;
                }else{
                    val = mu*servers;
                }
                dv += val;
                Q[0][rn-1] = val;
                Q[1][rn-1] = (i+1)-1;
                rn += 1;        
            }
            
            
            //add diagional element
            Q[0][rn-1] = -1*dv;
            Q[1][rn-1] = i+1;
            rn += 1;
            
            Q[2][i+1] = rn;
            
            //update the uniformization rate
            if(dv>gamma){
                gamma = dv;
            }
        }
        
        
        //System.out.println("Transition rate matrix completed.");
        
    }
    
  private void stochasticmatrix(){
        //converts Q to the stochastic matrix
      
        int rstart;
        int rstop;
        
        int i; int j; int r; int c;
        for(i=0; i<Ns; i++){
            r = i+1;
            rstart = (int) Q[2][i]; rstop = (int) Q[2][i+1];
            for(j=(rstart-1); j<(rstop-1); j++){
                c = (int) Q[1][j];
                if (c==r){
                    Q[0][j] = Q[0][j]*(1/gamma)+1;
                }else{
                    Q[0][j] = Q[0][j]*(1/gamma);
                }               
            }
            
        }
        
        //System.out.println("Converted to transition probability matrix.");
        
   }  
    
    
  private double numbiter(double t, double epsilon){
      //number of iterations in uniformization.
      //takes the time t, the tolerance epsilon, and the
      //uniformization rate gamma.
      
        double sigma = 1; double si = 1; double K = 0; 
        double tol = (1-epsilon)*Math.exp(gamma*t);
        while (sigma<tol){           
            si = si*((gamma*t)/(K+1));
            sigma = sigma + si;
            K += 1;
        }
            
        return K;
    }  
    
  
  public void uniformization(double t){
        //apply uniformization (also denoted randomization) for amount of time t
        
        double eps = 0.000001; //tolerance
        int K = (int) numbiter(t,eps);
      
        double[] y;
        double[] A = new double[Q[0].length];
        double[] z = new double[Ns];
        double gammat = gamma*t;
        
        int start;
        int stop;
        
        //initialize
        pi = copydoublevec(pi0);
        y = copydoublevec(pi0);
        
        //iterate
        int k; int i; int j; int c; int w;
        for(k=1; k<=K; k++){
 
            //regular scalar multiplication
            for(w=0; w<Q[0].length; w++){
                A[w] = Q[0][w]*(gammat/k);
            }
            
            //vector matrix multiplication
            z = fillout(z,0);
            for(i=0; i<Ns; i++){
                start = (int) Q[2][i];
                stop = (int) Q[2][i+1];
                for(j=(start-1); j<(stop-1); j++){
                    c = (int) Q[1][j];
                    z[c-1] += A[j]*y[i];
                }
            }
            y = copydoublevec(z);
            
            //pi + y
            pi = addvectors(pi,y);      
        }
        
        //finalize
        double constant = Math.exp(-1*gammat);
        pi = scalarmult(pi,constant);
        
        
    }
  
  public double blockingProbability(){
      
      double pB = pi[pi.length-1];
      
      return pB;
  }
  
  public double[] getStateDistribution(){
      
      int i;
      double[] stateDistOut = new double[pi.length];
      for(i=0; i<pi.length; i++){
          stateDistOut[i] = pi[i];
      }
      
      return stateDistOut;
  }
  
  
  public double expectedValue(){
        //calculates the expected value of the distribution dist
        double expVal = 0;
        
        int i;
        for(i=0; i<pi.length; i++){
            expVal += i * pi[i];
        }
        
        return expVal;
    }
    
    
 private static double[] copydoublevec(double[] x){
        
        int l = x.length;
        double[] y = new double[l];
        
        int i;
        for(i=0; i<l; i++){
            y[i] = x[i];
        }
        
        return y;
    }

private static double[] scalarmult(double[] x, double c){
        
        int l = x.length;
        int i;
        for(i=0; i<l; i++){
            x[i] = x[i]*c;
        }
        
        return x;
    }

private static double[] fillout(double[] v, double k){
        
        int l = v.length;
        int i;
        for(i=0; i<l; i++){
            v[i] = k;
        }
        
        return v;
    }

private static double[] addvectors(double[] x, double[] y){
        
        int l = x.length;
        double[] z = new double[l];
        int i;
        for(i=0; i<l; i++){
            z[i] = x[i]+y[i];
        }    
        return z;
    }
  
  
    
    
}
