FROM wpayne/openbabel:latest
MAINTAINER sanjaym@iiitd.ac.in 

# Install system dependencies
RUN apt-get update
RUN apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    python3.9 \
    python3-pip

RUN ln -s /usr/bin/python3.9 /usr/local/bin/python 
RUN ln -s /usr/bin/python3.9 /usr/local/bin/python3 
RUN rm -rf /var/lib/apt/lists/*


# Download and install Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh 
RUN /bin/bash /tmp/miniconda.sh -b -p /opt/conda 
RUN rm /tmp/miniconda.sh
RUN echo "export PATH=/opt/conda/bin:\$PATH" > /etc/profile.d/conda.sh 
RUN chmod +x /etc/profile.d/conda.sh



COPY Gcoupler.yml /tmp/Gcoupler.yml
RUN /bin/bash -c "source /etc/profile.d/conda.sh && conda env create -f /tmp/Gcoupler.yml && conda clean -afy"


# Set environment variables
ENV PATH="/opt/conda/envs/Gcoupler/bin:${PATH}"
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8


#RUN mkdir OB
RUN mkdir LG
#RUN mkdir cmake
RUN mkdir CDir

COPY LigBuilderV3.0.tar.gz /LG/
#COPY openbabel-2.4.1.tar.gz /OB/

# Install dependencies
RUN apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    gcc \
    make 
RUN rm -rf /var/lib/apt/lists/*

#WORKDIR /cmake/
#RUN wget https://cmake.org/files/v3.4/cmake-3.4.1-Linux-x86_64.tar.gz
#RUN tar -xf cmake-3.4.1-Linux-x86_64.tar.gz


RUN tar -xvf /LG/LigBuilderV3.0.tar.gz -C /LG/
WORKDIR /LG/LigBuilderV3.0/
RUN chmod 775 configure
RUN ./configure
WORKDIR /LG/LigBuilderV3.0/build/
RUN chmod +x build
RUN chmod +x cavity

WORKDIR /root
# Activate the environment
CMD ["conda", "init", "bash"]
CMD ["source", ".bashrc"]

WORKDIR /root
RUN echo "ulimit -s unlimited" >> .bashrc
CMD ["source", ".bashrc"]


CMD ["conda", "activate", "Gcoupler"]

RUN pip install py3Dmol
RUN pip install jupyter
RUN pip install -i https://test.pypi.org/simple/ Gcoupler


WORKDIR CDir
