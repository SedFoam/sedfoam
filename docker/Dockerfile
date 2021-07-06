FROM opencfd/openfoam2106-default
LABEL maintainer "cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10

RUN apt-get update \
  && apt-get install --no-install-recommends -y \
  git python3 python3-pip unzip python \
  mercurial libreadline-dev \
  texlive dvipng python3-tk \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
RUN pip3 install fluidfoam

USER 1000

WORKDIR /home/sudofoam
RUN git clone --recurse-submodules https://github.com/sedfoam/sedfoam.git

WORKDIR /home/sudofoam
RUN /bin/bash -c "hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -r vcompile/p2106 && cp sedfoam/docker/Allwmakeswak swak4Foam/Libraries/Allwmake"

WORKDIR /home/sudofoam/swak4Foam
RUN /bin/bash -c "shopt -s expand_aliases && ./maintainanceScripts/compileRequirements.sh"
WORKDIR /home/sudofoam/swak4Foam
RUN /bin/bash -c "shopt -s expand_aliases && source /usr/lib/openfoam/openfoam2106/etc/bashrc && ./Allwmake"

WORKDIR /home/sudofoam
RUN rm -rf swak4Foam

WORKDIR /home/sudofoam/sedfoam
RUN /bin/bash -c "shopt -s expand_aliases && source /usr/lib/openfoam/openfoam2106/etc/bashrc && ./Allwmake && ./Allwclean"
#RUN find /home/ofuser/sedfoam -name "*.o" -exec rm -f {} \;

WORKDIR /home/sudofoam
ENV HOME /home/sudofoam
ENTRYPOINT ["/bin/bash", "-c"]
