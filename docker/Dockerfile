FROM opencfd/openfoam-default:2112
LABEL maintainer "cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10

RUN update-ca-certificates && apt-get update \
  && apt-get install --no-install-recommends -y \
  git python3 python3-pip unzip python-is-python3 ipython3 \
  mercurial libreadline-dev vim nano emacs neovim \
  texlive dvipng python3-tk \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
RUN pip3 install --no-cache-dir fluidfoam pandas matplotlib odfpy

WORKDIR /root/
RUN git clone --branch develop --recurse-submodules https://github.com/sedfoam/sedfoam.git

WORKDIR /root
RUN /bin/bash -c "hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -b feature/port/p2112 && cp sedfoam/docker/Allwmakeswak swak4Foam/Libraries/Allwmake"

WORKDIR /root/swak4Foam
RUN /bin/bash -c "shopt -s expand_aliases && ./maintainanceScripts/compileRequirements.sh"

WORKDIR /root/swak4Foam
RUN /bin/bash -c 'shopt -s expand_aliases && source /openfoam/bash.rc && export FOAM_USER_LIBBIN=$FOAM_SITE_LIBBIN && export FOAM_USER_APPBIN=$FOAM_SITE_APPBIN && ./Allwmake'

WORKDIR /root/sedfoam
RUN /bin/bash -c 'shopt -s expand_aliases && source /openfoam/bash.rc && export FOAM_USER_LIBBIN=$FOAM_SITE_LIBBIN && ./Allwmake -prefix=group'

WORKDIR /home/sudofoam
RUN rm -rf swak4Foam sedfoam

USER 1000

WORKDIR /home/sudofoam
RUN git clone --branch develop --recurse-submodules https://github.com/sedfoam/sedfoam.git

WORKDIR /home/sudofoam
ENV HOME /home/sudofoam
ENTRYPOINT ["/bin/bash", "-c"]
