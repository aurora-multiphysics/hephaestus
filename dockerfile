# sudo docker build -t hephaestus .

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update
RUN apt-get install -y apt-utils
RUN apt-get install -y dialog
RUN apt-get install -y software-properties-common
RUN apt-get install -y git
RUN apt-get install -y make
RUN apt-get install -y cmake
RUN apt-get install -y libboost-all-dev
RUN apt-get install -y python3
RUN apt-get install -y libarmadillo-dev
RUN apt-get install -y clang-format-8
RUN apt-get install -y cppcheck
RUN apt-get install -y similarity-tester
RUN apt-get install -y flawfinder
RUN apt-get install -y doxygen
RUN apt-get install -y graphviz

RUN git clone https://github.com/terryyin/lizard.git /home/lizard
RUN cp /home/lizard/lizard.py /usr/local/bin
ENV PYTHONPATH /home/lizard

WORKDIR /home

CMD ["bash"]

