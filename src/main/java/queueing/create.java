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

/**
 * Generate queueing network transition rate matrix.
 * 
 * @author Anders Reenberg Andersen, andersra@live.dk
 */
public class create {
    
    double[][] Q; //transition rate matrix (in Harwell-Boeing format)
    
    double[][] rates; //adjacency matrix with rates (ser. rate only per server)
    double[][] Adj; //input adjacency (or routing) matrix
    double[] lambda; //arrival rate for each source node
    double[] mu; //service rate at each queue node 
    int[] servers; //number of servers at each queue node
    int[] cap; //capacity at each queue node
    double maxDiag; //maximum absolute diagonal element in Q
    boolean rejectWhenFull; //Can customers be rejected when downstream nodes are overcrowded
    
    
    int srcNodes; //number of source nodes (must correspond to lambda.length)
    int queueNodes; //number of queue nodes (must correspond to mu.length)
    
    public create(double[][] Adj, double[] lambda, double[] mu, 
            int[] servers, int[] cap, boolean rejectWhenFull){
        
        //required structure of Adj:
        //1. source nodes
        //2. queue nodes
        //3. (the) sink node
        
        //required structure of lambda, mu, c and cap:
        //values must occur in the same order as the nodes occur in Adj
        
        this.lambda = lambda;
        this.mu = mu;
        this.servers = servers;
        this.cap = cap;
        this.Adj = Adj;
        this.rejectWhenFull = rejectWhenFull;
        
        checkAdjancency();
        transitionRateMatrix();
        
    }
    
    private void checkAdjancency(){ //get network characteristics from the adjacency matrix 
       
        rates = new double[Adj.length][Adj[0].length];
        srcNodes = 0;
        queueNodes = 0;
        
        int n; int i;
        double sumRow; double sumCol;
        
        //determine node types
        for (n=0; n<(Adj.length-1); n++){
            sumRow = 0;
            for (i=0; i<Adj[n].length; i++){
                sumRow += Adj[n][i];
            }
            sumCol = 0;
            for (i=0; i<Adj.length; i++){
                sumCol += Adj[i][n];
            }
            
            if (sumRow>0&&sumCol==0){ //source node
                for (i=0; i<Adj[n].length; i++){
                    rates[n][i] += Adj[n][i]*lambda[srcNodes]; 
                }
                srcNodes++;
            }else if (sumRow>0&&sumCol>0){ //queue node
                for (i=0; i<Adj[n].length; i++){
                    rates[n][i] += Adj[n][i]*mu[queueNodes];
                }
                queueNodes++;    
            }else if (sumRow==0&&sumCol==0){
                System.err.println("Redundant node. Node neither sends nor receives.");
            }
        }
            
        
    }
    
    private void transitionRateMatrix(){
        
        //size of statespace
        int Ns = stateSpaceSize();
        //number of non-zero elements
        int nz = numberNonZero(Ns);
//        System.out.println("Number of states: " + Ns);
//        System.out.println("Non-zero elements: " + nz);
        
        Q = new double[3][];
        Q[0] = new double[nz]; //aa: Values
        Q[1] = new double[nz]; //ja: Column-index of values
        Q[2] = new double[Ns+1]; //ia: Range in aa and ja (value for each state + 1)
        
        //----------------------------------------
        
        int r; //row
        int c; //column
        int q; int q1; //queue
        double val = 0; //rate
        double diag; //diagonal
        int i = 0; //non-zero elements added
        int[] s = new int[queueNodes]; //state
        s[queueNodes-1] = -1;
        
        boolean sourced;
        boolean reject;
        maxDiag = 0;
        
        int[] qProd = new int[queueNodes];
        int prod = 1;
        for (q=(queueNodes-1); q>=0; q--){
            qProd[q] = prod;
            prod *= (cap[q]+1);
        }
        
        for(r=1; r<=Ns; r++){
            
            //get next state
            s = getNextState(s);
            
            Q[2][(r-1)] = i+1;
            
            diag = 0;
            for (q=0; q<queueNodes; q++){
              //sourcing
              sourced = false;
              for (q1=0; q1<srcNodes; q1++){
                if (Adj[q1][srcNodes+q]>0&&s[q]<cap[q]){
                    c = r + qProd[q];
                    val = rates[q1][srcNodes+q];
                    diag -= val;
                    Q[0][i] += val;
                    Q[1][i] = c;
                    sourced = true;
                    
                }
              }
              if (sourced){
                //System.out.println("Source " + Q[0][i]);
                i++;
              }
              
              
              
              //transfer and remove
              if (s[q]>0){  
                reject = false;
                for (q1=srcNodes; q1<Adj[0].length; q1++){ //receiving node
                        if ((q1-srcNodes)!=q&&Adj[srcNodes+q][q1]>0){ //basic demands
                            if (q1<(Adj[0].length-1)&&s[q1-srcNodes]<cap[q1-srcNodes]){ //transfer
                                c = r - qProd[q] + qProd[q1-srcNodes];
                                val = rates[srcNodes+q][q1]*min(s[q],servers[q]);
                                diag -= val;
                                Q[0][i] = val;
                                Q[1][i] = c;
                                //System.out.println("Transfer " + Q[0][i]);
                                i++;
                               
                            }else if(q1<(Adj[0].length-1)&&s[q1-srcNodes]==cap[q1-srcNodes]){
                                reject = true;
                            }else if (q1==(Adj[0].length-1)){ //remove
                                c = r - qProd[q];
                                val = rates[srcNodes+q][q1]*min(s[q],servers[q]);
                                diag -= val;
                                Q[0][i] = val;
                                Q[1][i] = c;
                                
                                if (rejectWhenFull==false||reject==false){
                                    //System.out.println("Remove " + Q[0][i]);
                                    i++;
                                }
                            
                                
                            }
                        }
                }
              
                //transfer-reject
                if (rejectWhenFull&&reject){
                    for (q1=srcNodes; q1<(Adj[0].length-1); q1++){
                        if ((q1-srcNodes)!=q&&Adj[srcNodes+q][q1]>0&&
                                s[q1-srcNodes]==cap[q1-srcNodes]){
                            
                                c = r - qProd[q];
                                val = rates[srcNodes+q][q1]*min(s[q],servers[q]);
                                diag -= val;
                                Q[0][i] += val;
                                Q[1][i] = c;
                                reject = true;
                        }
                    }
                    //System.out.println("Reject " + Q[0][i]);
                    i++;
                    
                }
              }
              
              
            }
            //diagonal
            Q[0][i] = diag;
            Q[1][i] = r;
            i++;
            
            //update maxDiag
            if (diag<maxDiag){
                maxDiag = diag;
            }
            
        }
        Q[2][Ns] = i+1;
        
        maxDiag = -1*maxDiag;
        
    }
    
    
    private int stateSpaceSize(){
        int Ns = 1; int i;
        for (i=0; i<cap.length; i++){
            Ns *= (cap[i]+1);
        }
        
        return Ns;
    }

    private int numberNonZero(int Ns){
        int nz = 0;
        
        //non-diagonal elements
        int i; int j; int k;
        int q; int l; int tr; int rm;
        for (i=0; i<(Adj.length-1); i++){
            q = 0; //receiving queue
            
            tr = 0; //number of receiving queues (transfer)
            rm = 0; //remove
            for (j=srcNodes; j<Adj[i].length; j++){
                if (Adj[i][j]>0){
                    if (i<srcNodes){ //i is a source node
                        l = cap[q];
                        for (k=0; k<queueNodes; k++){
                            if (k!=q){
                                l *= (cap[k]+1);
                            }
                        }
                        //System.out.println("Source: " + l);
                        nz += l;
                    }else{ //i is a queue node
                        if (j==(Adj[i].length-1)){ //discharge
                            rm++;
                            l = cap[(i-srcNodes)];
                            for (k=0; k<queueNodes; k++){
                                if (k!=(i-srcNodes)){
                                    l *= (cap[k]+1);
                                }
                            }
                            //System.out.println("Remove: " + l);
                            nz += l;
                        }else{ //transfer
                            tr++;
                            l = cap[(i-srcNodes)];
                            for (k=0; k<queueNodes; k++){
                                if (k!=(i-srcNodes)&&k!=q){
                                    l *= (cap[k]+1);
                                }
                            }
                            
                            if (rejectWhenFull){
                                l *= cap[q] + 1;
                            }else{
                                l *= cap[q];
                            }
                            
                            
                            //System.out.println("Transfer: " + l);
                            nz += l;
                            
                        }
                    }
                    
                }
                q++;
            }
            if (rejectWhenFull&&(rm+tr)>1){
                for (j=1; j<(rm+tr); j++){
                    nz -= cap[(i-srcNodes)]*j*binomialCoefficient((rm+tr),(j+1)); //subtract excess rejections
                }    
            }    
        }
        
        //add diagonal elements
        nz += Ns;
        
        return nz;
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
    
    
    private int min(int a, int b){
        if (a<b){
            return a;
        }else{
            return b;
        }
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
    
    
    public double[][] getTransitionRateMatrix(){
        //return the transition rate matrix
        
        return Q;
    }
    
    public double getUniformizationRate(){
        //returns the maximum absolute diagonal element
        //in the transition rate matrix.
        
        return maxDiag;
    }
    
            
    
}



