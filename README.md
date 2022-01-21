# hygroscopicityClosureTestbed

*Source code and data for "A Testbed to Evaluate Constraints on the Accuracy of Aerosol Mass Spectrometer Organic Aerosol Mass Fraction Measurements from Aerosol Hygroscopicity Data"*

# Data 

The ```data/``` directory contains files with all underlying data used in the manuscript.

# Building the testbed

The tesbed docker container is built via

```
docker build . -t mdpetters/hygroscopicitytestbed
```

A copy of the container is hosted on [DockerHub](https://hub.docker.com/repository/docker/mdpetters/hygroscopicitytestbed). To run the testbed install the Docker engine (https://docs.docker.com/engine/install/). Start the container:  

```bash
docker run -p 1234:1234 mdpetters/hygroscopicitytestbed
```

If the image is not present, docker will pull it. It takes a short time to initialize. To run the notebook go to

```
http://0.0.0.0:1234/open?path=webapp.jl
``` 
in your browser. 

# Manuscript code

The source Code to create the figures in the manuscript is provided in the ```manuscript/``` folder. The code was executed in Julia 1.7.1 using the provided environment. Note that the AIOMFAC executable must be present in the manuscript folder. You will need to compile the AIOMFAC model in the ```webapp/AIOMFAC.``` directory and copy the executable to ```manuscript/``` folder. For generic installation instructions on a linux computer please read the ```webapp/Dockerfile```. 


