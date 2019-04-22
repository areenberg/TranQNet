# Purpose
This library is specialized in evaluating any M/M/C/K queueing network that can be defined by an adjacency matrix, i.e. the model does not have to exist on *product form*. A weighted directed adjacency matrix is used by the class `create` to automatically construct the infinitesimal generator for the network at hand. The resulting object is then used as input in the class `evaluate` to retrieve the behavior of the network.

# Basic Overview

## Files

- `mc_math.jar`: The library file.

## Classes

- `create`: Automatically constructs the infinitesimal generator (the transition rate matrix) and various other parameters.

- `evaluate`: Uses an object defined with `create` to evaluate the queueing network.

# Getting Started
Firstly, download and add `mc_math.jar` to your Java-project.

Now import the `queueing` classes.

```
import queueing.*;
```

Define the weighted directed adjacency matrix using the structure: (1) Sources nodes, (2) queueing nodes, and (3) sink node. In this example, we have three queueing nodes. A single source node feeds all arriving customers into the first queueing node. The flow then splits into three parts sending 45% of the customers to queue 2, 30% to queue 3, and 25% to the sink node. Queue 2 and 3 send all customers to the sink after their service has been completed. 

```
double[][] A = {{0,1,0,0,0},
                {0,0,0.45,0.30,0.25},
                {0,0,0,0,1},
                {0,0,0,0,1},
                {0,0,0,0,0}};
```


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
