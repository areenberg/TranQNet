# Introduction
The purpose of this library is to provide a method for evaluating both transient and steady-state M/M/C/K single queues and queueing networks. The library has been designed to make this process as smooth as possible by restricting the input to the adjacency matrix and the basic parameters that characterize the system.

Using the library involves the following two steps:

**Step 1:** Construct the infinitesimal generator for the network at hand. This is handled automatically by the `create` class.

**Step 2:** Use the object created in the aforementioned step to retrieve the behavior of the network. This is handled by the `evaluate` class.

# Basic Overview

## Files

- `mc_math-1.2.jar`: The library file.

## Classes

- `create`: Automatically constructs the infinitesimal generator (the transition rate matrix) and various other parameters.

- `evaluate`: Uses an object defined with `create` to evaluate the queueing network.

# Getting Started
Firstly, download and add `mc_math.jar` to your Java-project.

Now import the `queueing` classes.

```java
import queueing.*;
```

Define the weighted directed adjacency matrix using the structure: (1) Sources nodes, (2) queueing nodes, and (3) sink node. In this example, we have three queueing nodes. A single source node feeds all arriving customers into the first queueing node. The flow then splits into three parts sending 45% of the customers to queue 2, 30% to queue 3, and 25% to the sink node. Queue 2 and 3 send all customers to the sink after their service has been completed.

```java
double[][] A = {{0,1,0,0,0},
                {0,0,0.45,0.30,0.25},
                {0,0,0,0,1},
                {0,0,0,0,1},
                {0,0,0,0,0}};
```

Define the remaining characteristics of the system, i.e. the arrival rate (`lambda`), service rates (`mu`), number of servers (`c`), capacity (`cap`), and how much of the capacity is occupied at time=0 (`occupiedCap`). If customers should be rejected when downstream queues are full, set `rejectWhenFull = true`; otherwise `rejectWhenFull = false`.  

```java
double[] lambda = {2}; double[] mu = {1.5,4,2.5}; int[] c = {2,1,2}; int[] cap = {20,20,20};
int[] occupiedCap = {0,0,0}; boolean rejectWhenFull = false;
```

Create the model.

```java
create network = new create(A,lambda,mu,c,cap,rejectWhenFull);
```

Prepare the evaluation calculations by plugging the `network` object into `evaluate`.

```java
evaluate system = new evaluate(occupiedCap,network);
```

Evaluate the behavior of the system at time=5 with a precision of 1x10^-9.

```java
system.uniformization(5,1e-9);
```

Evaluate the steady-state behavior of the system (i.e. at time=Inf) with a precision of 1x10^-6.

```java
system.gauss_seidel(1e-6);
```

Get the marginal state distribution for each of the network queues.

```java
double[][] dist = system.getMarginalDistributions();
```

Get the expected number of customers within each queue node.

```java
double[] expValue = system.expectedCustomers();
```

Get the probability of waiting on arrival at each queue node.

```java
double[] waitProb = system.waitingProbability();
```

# Citing this library

[![DOI](https://zenodo.org/badge/175442551.svg)](https://zenodo.org/badge/latestdoi/175442551)


Anders Reenberg Andersen. (2021). areenberg/mc_math: Converted to Maven project (v1.2). Zenodo. https://doi.org/10.5281/zenodo.5650104

# License
Copyright 2021 Anders Reenberg Andersen, PhD

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
