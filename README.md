# sphNG - dusty disc simulation

To run this first simulation you will need to checkout the branch named **dusty_sims**. This is a pretty complex simulation for the first for someone to run but the most expensive part is what we are likely to attempt to offload to GPU first.

Compile the code with,
> make openmp=yes gradhrkrt

Copy the executable **sph_tree_rk_gradh_RT** and the **inspho** file from the main sphNH folder to your working directory. This folder will also require the BonnorEbert.txt and setup.txt files, find them in the attched tar.

From the tar, you will also need to extract the folder called **Tables** to a location of your choice and set the environment variable **SPH_HOME** as the path to that folder.

To setup the initial conditions,
> ./sph_tree_rk_gradh_RT initial ifile < setup.txt out_setup.txt

Then we can start a simulation,
> ./sph_tree_rk_gradh_RT evolution ifile > stdout.txt