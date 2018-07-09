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

WORKDIR "/pycrysfml"

RUN ./build.sh

CMD /bin/bash
