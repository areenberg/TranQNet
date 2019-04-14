# Purpose
To easily and efficiently evaluate queueing networks based on numerical evaluations of their underlying Markov chain.

**Note:** *This branch is currently under heavy development, and the following information might not correspond to its current content.*

# Content
Currently this repository consists of a single Java-library (`mc_math.jar`) containing a single class. This class, denoted `evaluate`, can be used to numerically evaluate the state probability distribution of a transient queueing system with finite capacity, i.e. a queueing system of type M/M/C/K. The parameters of the system (arrival rate, service rate, number of servers and capacity) are homogeneous. The model then takes either the initial state probability distribution or the initial number of occupied servers. Uniformization (also denoted randomization and Jensen's method) is used to evaluate the probability distribution at time `t`. 

# Files

- `mc_math.jar`: Library containing the class `evaluate`.

- `src/queueing`: Folder containing the source code for `mc_math.jar`.  

# Usage

## Class `evaluate`
Calculates the state distribution at time `t` for an M/M/C/K queueing system.

### Constructors

- `evaluate(int ch, double l, double m, int c, int cap)`: Number of occupied servers as input.

- `evaluate(double[] initDist, double l, double m, int c, int cap)`: Initial state distribution as input.  

- `ch`: Number of currently occupied servers.
- `l`: Arrival rate. 
- `m`: Service rate.
- `c`: Total number of servers in the system.
- `cap`: Capacity of the system.
- `initDist`: Initial state probability distribution.

### Methods

- `void uniformization(double t)`: Calculates the state distribution at time `t` using uniformization.
- `double[] getStateDistribution()`: Returns the resulting state distribution.
- `double expectedValue()`: Returns the expected state.
- `double blockingProbability()`: Returns the blocking probability, i.e. the probability of attaining the last state in the state distribution.

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
