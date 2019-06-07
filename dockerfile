FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update

RUN apt upgrade -y

RUN apt install -y build-essential \
	gfortran \
	python \
	python-dev \
	python-numpy \
	python-scipy \
	python-matplotlib \
	python-pip \
	subversion \
	git \
	swig

RUN git clone https://github.com/bumps/bumps.git

RUN pip install ./bumps

RUN git clone https://github.com/scattering/pycrysfml.git
RUN git checkout origin/python3

WORKDIR "/pycrysfml"


RUN apt install -y python3-pip \
	build-essential libssl-dev libffi-dev python3-dev
	2to3

RUN ./build.sh

RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install bumps

CMD /bin/bash
