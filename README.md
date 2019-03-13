# Purpose
The easily and efficiently evaluate queueing models based on numerical evaluations of their underlying Markov chains.

# Content
Currently this repository consists of a single Java-library (`mc_math.jar`) containing a single class. This class, denoted `evaluate`, can be used to numerically evaluate the state probability distribution of a transient queueing system with finite capacity, i.e. a queueing system of type M/M/C/K. The parameters of the system (arrival rate, service rate, number of servers and capacity) are homogeneous. The model then takes either the initial state probability distribution or the initial number of non-idle servers. Uniformization (also denoted randomization and Jensen's method) is used to evaluate the probability distribution at time `t`. 

# Files

- `mc_math.jar`: Library which contains the `evaluate` class.

- `src/queueing`: Folder containing the source code for `mc_math.jar`.  

# Usage

## Class `evaluate`

### Constructors

- `evaluate(int ch, double l, double m, int c, int cap)`: Number of non-idle servers as input.

- `evaluate(double[] initDist, double l, double m, int c, int cap)`: Initial distribution as input.  

- `ch`: Number of servers currently in use.
- `l`: Arrival rate. 
- `m`: Service rate.
- `c`: Total number of servers in the system (idle and non-idle).
- `cap`: Capacity of the system.
- `initDist`: Initial state probability distribution.

### Methods

- `uniformization(double t)`: Apply uniformization at time `t`.
- `getStateDistribution()`: Returns the resulting probability state distribution.
- `expectedValue()`: Returns the expected state.
- `blockingProbability()`: Returns the blocking probability, i.e. the probability of attaining the last state in the state distribution.

# License
Copyright 2019 Anders Reenberg Andersen, PhD

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
