/*
 *Copyright 2021 Anders Reenberg Andersen, PhD
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
  int[] servers; //the number of servers
  double[][] Q;
  double[] pi0; //initial state probability distribution
  double[] pi; //resulting state probability distribution
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
      this.servers = network.servers;
      
      Ns = Q[2].length-1;
  }
  
  
  public evaluate(int[] occupiedCap, create network){
      //constructur for the custom queueing system with
      //current amount of occupied capacity as input.
      
      this.gamma = network.maxDiag;
      this.Q = network.Q;
      this.queueNodes = network.queueNodes;
      this.cap = network.cap;
      this.servers = network.servers;
      
      pi0 = new double[Q[2].length-1];
      int n = 1; int prod = 1;
      int i; 
      for (i=(occupiedCap.length-1); i>=0; i--){
          n += occupiedCap[i]*prod;
          prod *= (network.cap[i]+1);
      }
      pi0[n-1] = 1;
      
      Ns = Q[2].length-1; 
      
  }
  
  
  public evaluate(int occupiedCap, double lambda, double mu, int ser, int K){ 
      //constructor for the M/M/C/K model with occupied servers as input 
     
      
      this.currentlyUsed = occupiedCap;
      this.lambda = lambda;
      this.mu = mu;
      servers = new int[1];
      servers[0] = ser;
      
      Ns = K + 1;
      
      transitionratematrix();
      initialDistribution();
  }
  
  
  public evaluate(double[] initDist, double lambda, double mu, int ser, int K){ 
       //constructor for the M/M/C/K model with state distribution as input
      
      this.lambda = lambda;
      this.mu = mu;
      servers = new int[1];
      servers[0] = ser;
      
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
                if (v<=servers[0]){
                    val = mu*v;
                }else{
                    val = mu*servers[0];
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
    
    
  private double numbiter(double gammat, double epsilon){
      //number of iterations in uniformization.
      //takes the time t, the tolerance epsilon, and the
      //uniformization rate gamma.
      
        double sigma = 1; double si = 1; double K = 0; 
        double tol = (1-epsilon)*Math.exp(gammat);
        while (sigma<tol){           
            si = si*((gammat)/(K+1));
            sigma = sigma + si;
            K += 1;
        }
            
        return K;
    }  
    
  
  public void uniformization(double t, double eps){
        //apply uniformization (also denoted randomization) for amount of time t
        //and tolerance eps (e.g. t=1 and eps=1e-6).
        
        double[] P = stochasticmatrix();
        
        //calculate if gamma*t causes underflow
        double tUnderflow = 70.0/gamma;
        int steps = 1;
        double[] tvec = {t}; 
        if (t>tUnderflow){ //if gamma*t causes underflow, use uniformization in sequence
            steps = (int) Math.ceil(t/tUnderflow);
            tvec = new double[steps];
            for (int i=0; i<(steps-1); i++){
                tvec[i] = tUnderflow;
            }
            tvec[(steps-1)] = t - tUnderflow*(steps-1);
        }
            
        double[] y;
        double[] A = new double[P.length];
        double[] z = new double[Ns];
        double gammat;
        double constant;
        
        int start;
        int stop;
        int k; int K; int i; int j; int c; int w; int stp;
        
        //initialize
        pi = copydoublevec(pi0);
        y = copydoublevec(pi0);
        
        for (stp=0; stp<steps; stp++){ 
        
            //get number of iterations
            gammat = gamma*tvec[stp];
            K = (int) numbiter(gammat,eps);
        
            //iterate
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
            
                pi = addvectors(pi,y);      
            }
        
            //finalize
            constant = Math.exp(-1*gammat);
            pi = scalarmult(pi,constant);
        
            if (stp<(steps-1)){
                y = copydoublevec(pi);
            }
        
        }
        
  }
  
  public void gauss_seidel(double tol){
        //calculates the steady-state solution using Gauss-Seidel
        
        int z;
        int qw;
        int sw;
        int lw;
        
        //transpose Q
        double[][] A = transpose();
        
        //scale so aii = 1 for all i 
        A = scale(A);
        
        //start algorithm      
        double startTime = System.currentTimeMillis();
        
        double sum;
        int jaj;
        int iter;
        double ss;
        
        double eps = 1;
        int m = 5;
        int k = 0;
        double[] pi_old;
        
        //initialize state distribution
//        pi = new double[Ns];
//        double ns = Ns;
//        for(z=0; z<Ns; z++){
//            pi[z] = 1.0/ns;
//        }
        uniformization(10,1e-3);
        
        
        while (eps>tol || k > 10000){
        pi_old = copydoublevec(pi);
            
        //run sequence
        for(iter=0; iter<m; iter++){
            k++;
        
            for(z=0; z<Ns; z++){
                sum = 0;
                sw = (int) A[2][z]-1;
                lw = (int) A[2][z+1]-1;
                for(qw= sw; qw<lw; qw++){
                    jaj = (int) A[1][qw];
                    sum = sum + A[0][qw]*pi[jaj-1];
                }
                pi[z] = pi[z] - sum;
                
                if (pi[z]<1e-14){
                    pi[z] = 0.0; //prevent underflow
                }
            }
            
        }
        
        
        //normalize                               
        ss = sumoverdouble(pi);
        for(z=0; z<Ns; z++){
            pi[z] = pi[z]/ss;
        }
       
        //check convergence
        eps = relativetol(pi_old, pi);
        m = adjustspace(m,eps);
        
        }
        
        
        double stopTime = System.currentTimeMillis();
        double elapsedTime = stopTime - startTime;
        //System.out.println("Solution converged in " + k + " iterations and " + (elapsedTime/1000) + " seconds");
     
//        double[] b = check(A,pi);
//        for(z=0; z<6; z++){
//            System.out.println("b: " + b[z]);
//        }
//        double sumb = sumoverdouble(b);
//        System.out.println("The sum over b is " + sumb);
        
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
  
  public double[] expectedCustomers(){
        //calculates the expected number of customers that are in the queueing system
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
 
  
  public double[] expectedQueueLength(){
        //calculates the expected length of the queue
        double[] expVal = new double[queueNodes];
        double[][] dist = getMarginalDistributions();
        
        int q; int i; int k;
        for (q=0; q<queueNodes; q++){
            if (servers[q]<cap[q]){
                k = 0;
                for(i=(servers[q]+1); i<dist[q].length; i++){
                    k++;
                    expVal[q] += k*dist[q][i];
                }
            }
        }
        
        return expVal;
    }
  
  public double[] waitingProbability(){
        //calculates the probability of waiting on arrival
        double[] waitProb = new double[queueNodes];
        double[][] dist = getMarginalDistributions();
        
        int q; int i; double sum;
        for (q=0; q<queueNodes; q++){
            sum = 0;
            for(i=0; i<servers[q]; i++){
                sum += dist[q][i];
            }
            waitProb[q] = 1 - sum;
        }
        
        return waitProb;
    }

  
  public double[] blockingProbability(){
        //calculates the probability of a customer being rejected on arrival
        double[] blockProb = new double[queueNodes];
        double[][] dist = getMarginalDistributions();
        
        int q;
        for (q=0; q<queueNodes; q++){
            blockProb[q] = dist[q][dist[q].length-1];
        }
        
        return blockProb;
    }
  
  public double[] emptyProbability(){
        //calculates the probability of a customer being rejected on arrival
        double[] emptyProb = new double[queueNodes];
        double[][] dist = getMarginalDistributions();
        
        int q;
        for (q=0; q<queueNodes; q++){
            emptyProb[q] = dist[q][0];
        }
        
        return emptyProb;
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
  

private double sumoverdouble(double[] M){
    double s = 0;
    int l = M.length;
    int i;
    for(i=0; i<l; i++){
        s = s + M[i];
    }
    
    return s;
}


private double relativetol(double[] x_old, double[] x_new){
    
    int l = x_old.length;
    
    double[] diff = new double[l];
    
    double largest = 0;
    int i;
    for(i=0; i<l; i++){
        diff[i] = Math.abs(x_new[i]-x_old[i])/Math.abs(x_old[i]);
        if(diff[i]>largest){
            largest = diff[i];
        }
    }
    
    return largest;
} 

private int adjustspace(int m_old, double eps){
         
         int m;
         double u = eps/m_old;
         
         m = (int) Math.round(45*Math.exp(-1*u))+5;
         
         return m;
}

private int binomialCoefficient(int n, int k){
        //From The flying keyboard
        // 2018 TheFlyingKeyboard and released under MIT License
        // theflyingkeyboard.net
        
        int top = 1;
        for(int i = n; i > k; --i){
            top *= i;
        }
        return top / factorial(n - k);
}
    
private int factorial(int n){
        //From The flying keyboard
        // 2018 TheFlyingKeyboard and released under MIT License
        // theflyingkeyboard.net
        
        int fact = 1;
        for(int i = 2; i <= n; i++){
            fact *= i;
        }
        return fact;
}

private double[][] transpose(){
    //transpose the transition rate matrix
    
    double[][] A = new double[3][];
    A[0] = new double[Q[0].length]; //aa
    A[1] = new double[Q[1].length]; //ja
    A[2] = new double[Q[2].length]; //ia
    
    //maximum number of transitions from s 
    int maxCol = queueNodes*2 + binomialCoefficient(queueNodes,2)*2 + 1; 
    
    double[][][] B = new double[Ns][2][maxCol]; //row,{value,column},nr
    int[] Bb = new int[Ns]; //number added to row in B
    
    int i; int j;
    int a; int b; int ri;
    //sort elements of Q into rows of A
    for (i=0; i<Ns; i++){
        a = (int) Q[2][i]; b = (int) Q[2][i+1];
        for (j=a; j<b; j++){
            ri = (int) Q[1][j-1] - 1; 
            B[ri][0][Bb[ri]] = Q[0][j-1];
            B[ri][1][Bb[ri]] = i + 1;
            Bb[ri]++;
        }
    }
    
    int sum = 0;
    for (i=0; i<Bb.length; i++){
        sum += Bb[i];
    }
    
    //put B into A
    int k;
    A[2][0] = 1;
    for (i=0; i<Ns; i++){
        A[2][i+1] = A[2][i] + Bb[i];
        a = (int) A[2][i]; b = (int) A[2][i+1];
        k = 0;
        for (j=a; j<b; j++){    
            A[0][j-1] = B[i][0][k];
            A[1][j-1] = B[i][1][k];
            k++;
        }
    }
    
    return A;
}


private double[][] scale(double[][] A){
    //scale matrix A so all diagonal elements equal 1
    
    double scale;
    int z; int sw; int lw;
    int qw; int slen = A[2].length-1;
    for(z=0; z<slen; z++){
        sw = (int) A[2][z]-1;
        lw = (int) A[2][z+1]-1;
        scale = 0;
        for(qw=sw; qw<lw; qw++){
            if(A[0][qw]<scale){
                scale = A[0][qw];
            }    
        }
        scale = 1.0/scale;
        for(qw= sw; qw<lw; qw++){
            A[0][qw] = A[0][qw]*scale;
        }
    }
    
    return A;
}


private double[] check(double[][] A, double[] x){
        double[] z = new double[Ns];
        
        double sum;
        int initial;
        int last;
        int i;
        int j;
        int jaj;
        
        for(i=0; i<Ns; i++){
            sum = 0;
            initial = (int) A[2][i];
            last = (int) A[2][i+1]-1;
            for(j= initial; j<=last; j++){
                jaj = (int) A[1][j-1];
                sum = sum + A[0][j-1]*x[jaj-1];
            }
            z[i] = sum;
        }
        
        
        return z;
    }



}
