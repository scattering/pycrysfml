FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y build-essential \
	gfortran \
	python3 \
	python3-dev \
	python3-numpy \
	python3-scipy \
	python3-matplotlib \
	python3-pip \
	subversion \
	git \
	swig

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1
RUN pip3 install git+https://github.com/bumps/bumps.git

RUN git clone https://github.com/scattering/pycrysfml.git

WORKDIR "/pycrysfml"
RUN ./build.sh


CMD /bin/bash
