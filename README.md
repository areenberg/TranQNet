# Purpose
Evaluation of queueing models based on numerical evaluations of Markov chains.

# Content
Currently this repository consists of a single Java-library containing a single class. This class, denoted `evaluate`, can be used to numerically evaluate the state probability distribution of a transient queueing system with finite capacity, i.e. a queueing system of type M/M/C/K. The parameters of the system (arrival rate, service rate, number of servers and capacity) are homogeneous. The model then takes either the initial state probability distribution or the initial number of non-idle servers. Uniformization (also denoted randomization and Jensen's method) is used to evaluate the probability distribution after time t. 

# Files and Folders

- `mc_math.jar`: Library which contains the `evaluate` class.

- `src/queueing`: Contains the source code associated with `mc_math.jar`.  

# Usage




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
