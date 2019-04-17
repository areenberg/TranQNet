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
  int queueNodes; //number of queue nodes in the system
  int[] cap; //capacity at each queue node
  
  public evaluate(double[] initDist, create network){
      //constructur for the custom queueing system with 
      //state distribution as input.
      
      this.pi0 = initDist;
      this.gamma = network.maxDiag;
      this.Q = network.Q;
      this.queueNodes = network.queueNodes;
      this.cap = network.cap;
      
      Ns = Q[2].length-1;
  }
  
  
  public evaluate(int[] occupiedServers, create network){
      //constructur for the custom queueing system with
      //current number of occupied servers as input.
      
      this.gamma = network.maxDiag;
      this.Q = network.Q;
      this.queueNodes = network.queueNodes;
      this.cap = network.cap;
      
      pi0 = new double[Q[2].length-1];
      int n = 1; int prod = 1;
      int i; 
      for (i=(occupiedServers.length-1); i>=0; i--){
          n += occupiedServers[i]*prod;
          prod *= (network.cap[i]+1);
      }
      pi0[n-1] = 1;
      
      Ns = Q[2].length-1; 
      
  }
  
  
  public evaluate(int occupiedServers, double lambda, double mu, int servers, int K){ 
      //constructor for the M/M/C/K model with occupied servers as input 
     
      
      this.currentlyUsed = occupiedServers;
      this.lambda = lambda;
      this.mu = mu;
      this.servers = servers;
      
      Ns = K + 1;
      
      transitionratematrix();
      initialDistribution();
  }
  
  
  public evaluate(double[] initDist, double lambda, double mu, int servers, int K){ 
       //constructor for the M/M/C/K model with state distribution as input
      
      this.lambda = lambda;
      this.mu = mu;
      this.servers = servers;
      
      Ns = K + 1;
      
      transitionratematrix();
      
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
        //generates the transition rate matrix for an M/M/C/K queueing system in Harwell-Boeing format
        
      
        int nz = Ns+(Ns-1)*2; //non-zero elements (capacity x 2 + diagonal)
        
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
                rn++;        
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
                rn++;        
            }
            
            
            //add diagional element
            Q[0][rn-1] = -1*dv;
            Q[1][rn-1] = i+1;
            rn++;
            
            Q[2][i+1] = rn;
            
            //update the uniformization rate
            if(dv>gamma){
                gamma = dv;
            }
        }
        
        
    }
    
  private double[] stochasticmatrix(){
        //converts Q to the stochastic matrix
      
        double[] P = new double[Q[0].length];
        
        int rstart;
        int rstop;
        
        int i; int j; int r; int c;
        for(i=0; i<Ns; i++){
            r = i+1;
            rstart = (int) Q[2][i]; rstop = (int) Q[2][i+1];
            for(j=(rstart-1); j<(rstop-1); j++){
                c = (int) Q[1][j];
                if (c==r){
                    P[j] = Q[0][j]*(1/gamma)+1;
                }else{
                    P[j] = Q[0][j]*(1/gamma);
                }               
            }
            
        }
        
        return P;
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
    
  
  public void uniformization(double t, double eps){
        //apply uniformization (also denoted randomization) for amount of time t
        //and tolerance eps (e.g. t=1 and eps=0.000001).
        
        double[] P = stochasticmatrix();
        
        int K = (int) numbiter(t,eps);
      
        double[] y;
        double[] A = new double[P.length];
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
            for(w=0; w<P.length; w++){
                A[w] = P[w]*(gammat/k);
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
  
  
  public double[] getStateDistribution(){
      //return the state distribution
      
      return pi;
  }
  
  public double[][] getMarginalDistributions(){
      //calculates and returns the marginal probability
      //distributions for each queue in the system.
      
      double dist[][] = null;
      int[] state;
      if (queueNodes>1){
          dist = new double[queueNodes][];
          int q; int s;
          for (q=0; q<queueNodes; q++){
              dist[q] = new double[cap[q]+1];
              state = new int[queueNodes];
              state[queueNodes-1] = -1;
              for (s=0; s<pi.length; s++){
                  state = getNextState(state);
                  dist[q][state[q]] += pi[s];
              }
          }
      }else{
          dist = new double[1][];
          dist[0] = pi;
      }
      
      return dist;
  }
  
  public double[] expectedValue(){
        //calculates the expected value for each queue in the system
        double[] expVal = new double[queueNodes];
        double[][] dist = getMarginalDistributions();
        
        int q; int i;
        for (q=0; q<queueNodes; q++){
            for(i=0; i<dist[q].length; i++){
                expVal[q] += i*dist[q][i];
            }
        }
        
        return expVal;
    }
    
    
 private double[] copydoublevec(double[] x){
        
        double[] y = new double[x.length];
        
        int i;
        for(i=0; i<x.length; i++){
            y[i] = x[i];
        }
        
        return y;
    }

private double[] scalarmult(double[] x, double c){
        
        int l = x.length;
        int i;
        for(i=0; i<l; i++){
            x[i] = x[i]*c;
        }
        
        return x;
    }

private double[] fillout(double[] v, double k){
        
        int l = v.length;
        int i;
        for(i=0; i<l; i++){
            v[i] = k;
        }
        
        return v;
    }

private double[] addvectors(double[] x, double[] y){
        
        double[] z = new double[x.length];
        int i;
        for(i=0; i<x.length; i++){
            z[i] = x[i]+y[i];
        }    
        return z;
    }


private int[] getNextState(int[] state){
        //returns the next state using the current state as input.
        
        int q = queueNodes-1; //move backwards from queue number
        while(state[q]==cap[q]){
            state[q] = 0;
            q--;
        }
        state[q]++;
        
        return state;
    }
  
  
    
    
}
