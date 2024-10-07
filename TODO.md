TODO
====

Or possible future upgrades. 


Parallelization using OpenMP. 

Now we have somehow homebrew parallelization using xargs and 
a program that does one piece at a tame. It is ok-ish, but could 
be improved with OpenMP. I have tried it at some point woth other 
program using capdDDEs, but I have failed. It seems capd
itself works ok with OpenMP, so it might be an issue with my code. 

The exemplary program is to try to implement for pieces is:

```
#include <iostream>
#include <omp.h>

int main() {
    // Set the number of threads
    int num_threads = 4;
    
    // Enable parallelism using OpenMP
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < 10; i++) {
    	 // do i-th piece here.
        int thread_id = omp_get_thread_num();
        std::cout << "Iteration " << i << " is being executed by thread " << thread_id << std::endl;
    }

    return 0;
}
```