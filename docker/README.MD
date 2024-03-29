# Docker container

<br>
<div align="center">
<img src="/Images/Docker_page.png"></div>
<br>

This is the repo of the official [Docker image](https://hub.docker.com/r/sanjayk741/gcoupler) for **Gcoupler**.


Instructions for installing and running docker on any PC can be found [here](https://docs.docker.com/engine/install/) 
1. [Windows](https://docs.docker.com/desktop/install/windows-install/)
2. [MacOS](https://docs.docker.com/desktop/install/mac-install/)
3. [Linux](https://docs.docker.com/desktop/install/linux-install/)

Pull the **Gcoupler** image from Docker Hub by running the following command in your terminal:
```
$ docker pull sanjayk741/gcoupler
```
Verify the new image has been created using the **docker images** command.
```
$ docker images
```
To access the terminal of a Docker image, you can use the **docker run** command with the **-it** option.
```
$ docker run -it <image-name> /bin/bash
```
Replace **<image-name\>** with the name or ID of the Docker image of **Gcoupler**.

Find the ID of the currently running container for **input** and **output**.
```
$ docker ps -a
```
To start the container again and access its terminal. 
```
$ docker start <container-ID>
$ docker exec -it <container-ID> bash
```

## Running **Gcoupler**
You can either run Gcoupler by python **command line** mode directly
```
$ python
>>> import Synthesizer
>>>
```

### OR 
By **Jupyter Notebook** in following 2 step

STEP-I

Run the following command to get the IP address of your running container
```
hostname -I
```
STEP-II

Launch Jupyter Notebook by running the following command:
```
$ jupyter-notebook --ip=<ip-address> --allow-root
```
Replace **<ip-address\>** with the IP address that you noted down in **STEP-I**

  This should start the Jupyter Notebook server and display a URL in the terminal that you can use to access the Jupyter Notebook interface.

### Default paths inside the container
```
>>> import Synthesizer as as
>>> sz.Set_paths(LigBuilder_path='/LG/LigBuilderV3.0/',libiomp5_file='/opt/conda/envs/Gcoupler/lib/',Output_dir='/root/CDir/')
```
The option **Output_dir** can be changed in case user decides to make other directories inside the container. 

## Input/Output

### Input 
Find the ID of the currently running container, just executed using the **docker ps -a** command.
```
$ docker ps -a
```
To write a file to the container, use the **docker cp** command to copy the file from the host to the container.
```
$ docker cp file container_id:/root/CDir/
```
This command will copy the **file** file (pdb/tsv) from the host's current directory to the **Gcoupler** container with ID **container_id** at the CDir/ directory inside the container.

### Output
Find the ID of the currently running container, just executed using the **docker ps -a** command.
```
$ docker ps -a
```
To write a file from the container, use the **docker cp** command to copy the file from the container to the host.
```
$ docker cp container_id:/root/CDir/file-name .
```
This command will copy the **file** file (csv/pdf)from the **Gcoupler** container with ID **container_id** under the CDir/ directory inside the container to the host's current directory.
