FROM ubuntu:16.04

MAINTAINER Konstantin Malanchev <malanchev@physics.msu.ru>

RUN apt-get update &&\
    apt-get install -y &&\
        gcc \
        g++ \
        cmake \
        libboost-all-dev \
        libjpeg-dev \
        libglew-dev \
        libglfw3-dev \
        libglm-dev \
        &&\

CP . /app
WORKDIR /app

RUN cmake --build cmake-build-debug --target main -- -j 4

CMD ["cmake-build-debug/main"]