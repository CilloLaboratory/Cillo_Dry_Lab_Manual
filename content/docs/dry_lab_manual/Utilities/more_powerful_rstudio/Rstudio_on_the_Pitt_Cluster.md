# Running RStudio on the Pitt HTC

Sometimes, we will need to use more resources for computation than are available on our laptops. When more CPUs or RAM are required, we can utilize the resources available to our lab through the [Center for Research Computing](https://crc.pitt.edu/resources/computing-hardware) or CRC. Specifically, our lab has an allocation on the Highthroughput Computing Cluster or HTC. More info on the hardware infrastructre at the CRC is available [here](https://crc.pitt.edu/resources/computing-hardware) for the interested reader.

While we have plenty of resources for analysis on the HTC, the resources are not unlimited and should only be used when more power is needed. As such, don't be afraid to use these resources but please use them judisciously. For example, do not run jobs for multiple days - even though it's possible, it's a terrible waste of resources!

# Goals of the process 
- Setup a reproducible Rstudio environment on the HTC using singularity
- Run a script to start an Rstudio instance on the HTC
- Login to the Rstudio instance on the HTC 
- Download the necessary packages in a local directory

We note Rstudio instances are available in a pre-built environment called OnDemand on the HTC. However, this has the downside of not being the same Rstudio environment that we will be using locally on our laptops and thus prevents us from customizing the environment like we can do with singularity. 

## Prerequsites 
For this tutorial, we will need...
    - Some familiarity with Linux command line operations, but we will cover what we need for the task at hand. Additional background and resources for basic familarity can be found [here](https://crc.pitt.edu/linux)
    - Some familiarity with SLURM, which can be found in [this blog](https://blog.ronin.cloud/slurm-intro/)
    - A user account on the Pitt HTC
    - A singularity image of rstudio (available in the Cillo Lab shared folder on OneDrive)
    - An sbatch script to start rstudio and specify resources (available in the Cell Lab shared folder on OneDrive)

## Log in to the HTC 
First, we will log in the the HTC. We do this by opening terminal on our Mac (or a terminal-like environment on Windows). If you are logging in off-campus, you will need to use a VPN client to do as as described [here](https://services.pitt.edu/TDClient/33/Portal/KB/ArticleDet?ID=293). Next, we will type the following: 

```tpl
ssh username@htc.crc.pitt.edu
```

Followed by inputting your Pitt password. Note that when you type your password, nothing will appear. This is expected for security reasons!

If you are becoming a pro and are logging in to the HTC a lot, you can set up [SSH keys](https://goteleport.com/blog/how-to-set-up-ssh-keys/) so you don't have to type your password every time.

## Setup a reproducible Rstudio environment on the HTC 
In a previous tutorial, we introduced Docker and explained how it is helpful for reproducibility. One downside of Docker is that it cannot be used in shared cluster computing environments because of security issues. Instead, we can use a similar platform called singularity. 

To start using singularity on the HTC, we need an image similar to what we have in Docker. It is beyond the scope of this tutorial to introduced how to create images for singularity, so instead we will rely on the pre-build image *cillo_lab_docker_v1.2.sif*. This image is available on the Cillo Lab shared folder on OneDrive.

Create a directory in your HOME directory (the directory where you start when you login) called "rstudio_singularity" and change into that directory with the following code:

```tpl
mkdir rstudio_singularity
cd rstudio_singularity
```

Upload a copy the *cillo_lab_docker_v1.2.sif* image to this directory. This can be done through [scp](https://www.linode.com/docs/guides/how-to-use-scp/) or potentially more easily through the OnDemand file system. 

Now, we are ready to create an Rstudio instances on the HTC.

## Start an Rstudio instance with an sbatch script
We can now start an Rstudio session with as many CPUs and gigabytes of RAM as well need for analysis. We will achieve this using an sbatch script that takes a few variables as input. This script is available [here](path/to/script).

We can use the following code in conjunction with the sbatch script mentioned above to start an Rstudio session with 2 CPUs for 1 hour:

```tpl
sbatch --cpus-per-task=2 --time=01:00:00 --job-name=rstudio --export=workdir=${PWD},r_lib_dir=${PWD}/Rlibs_arc85_Jan_2024 singularity_rstudio.sbatch
```

Executing this code should produce an output indicating that job is submitted. You can check the status with the following, where you can substitute your username for mine (e.g. replace arc85 with your Pitt username):

```tpl
squeue -u arc85
```

Finally, we can check the slurm output for some info about how to login to the running studio session via our web browser.


## Login to the HTC Rstudio from your web browser 

We will next open a port on our laptop to connect to the HTC. This is accomplished be copying the ssh code from the output of the slurm file into the Terminal on your laptop. Once you do this, it will pause that terminal window, indicating that the port is open. 

After the port is open, we will go to our web browser of choice (i.e. Google Chrome). We can then connect via the open port to the running Rstudio session on the HTC. Once you do this, you will be taken to a login screen for Rstudio. Input your username and password from the slurm output file on the HTC, and you will be logged in and ready to compute!

## Check your library path 

In rstudio, we can now check to be sure the packages we download will go to where we want them to. We can check this with:

```tpl
.libPath()
```

Which should point to the directory you specified upon startup. Importantly, this will allow the packages to persist from one session to another even through the RStudio session will be killed.

## Download packages

We can now download packages as usual from bioconductor and CRAN or wherever else we choose. These will be installed into the packages directory you specified and will be available for future use even when we stop this Rstudio server session.

## Proceed with analysis

Analysis can now proceed as if we were using our local version of Rstudio. The nice thing about singularity is that we can also write to the directory, so we can easily save files to the HTC for future use. We can even compile RMD files!

## Shutting down 

Whenever we are finished with our analysis on the HTC, we can shut everything down with these steps: 

- Power down Rstudio session in the browser
    - DO NOT save session .RData info when prompted
- Close the window
- CTL+C to kill the SSH 
- scancel (using the job number provided in the slurm outfile) to kill the resources
    - This is only required if there is still time remaining on your sbatch job