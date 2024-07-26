FROM ubuntu:jammy

RUN apt-get update 
RUN apt-get -y install apt-transport-https software-properties-common wget openssh-client
 
RUN wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc
RUN add-apt-repository http://dl.openfoam.org/ubuntu 
 
RUN wget -O - http://dl.openfoam.org/gpg.key | apt-key add - && add-apt-repository http://dl.openfoam.org/ubuntu

RUN apt-get update 
RUN apt-get install -y --no-install-recommends openfoam12 make 
RUN rm -rf /var/lib/apt/lists/*

RUN apt-get update 
RUN apt-get -y install bash-completion bash-builtins libnss-wrapper vim nano tree curl
RUN rm -rf /var/lib/apt/lists/*

RUN apt-get update 
RUN apt-get install -y emacs gedit gedit-plugins gnuplot gnuplot-x11 gnuplot-doc meld less
RUN apt-get install -y less unzip
RUN rm -rf /var/lib/apt/lists/*

RUN apt-get install -y make

RUN apt-get update && apt-get install -y python3 python3-pip python-is-python3
RUN pip3 install numpy shapely numpy-stl pandas scipy 

# add user and create group
RUN adduser --disabled-password openfoam

USER openfoam
WORKDIR /home/openfoam

RUN echo "USER=openfoam" >> .bashrc
RUN echo "source /opt/openfoam12/etc/bashrc" >> .bashrc
RUN /bin/bash -ic "source ~/.bashrc"
RUN curl -LOk https://github.com/demichie/OpenPDAC-12/archive/master.zip && unzip *.zip && rm master.zip
RUN /bin/bash -ic "/home/openfoam/OpenPDAC-12-main/applications/OpenPDAC/Allwmake"
RUN mkdir -p /home/openfoam/OpenFOAM/openfoam-12

