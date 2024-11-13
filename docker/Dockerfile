FROM opencfd/openfoam-default:2406
LABEL maintainer="cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10
ARG DEBIAN_FRONTEND=noninteractive

RUN update-ca-certificates && apt-get update \
  && apt-get install --no-install-recommends -y \
  python3-pip python3-dev python-is-python3 ipython3 \
  unzip git mercurial libreadline-dev vim nano emacs python3-venv \
  #python3-tk cm-super meld \
  #texlive dvipng texlive-latex-extra python3-netcdf4 python3-pandas \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /home/sudofoam/pyenv

RUN /bin/bash -c 'shopt -s expand_aliases && source /home/sudofoam/pyenv/bin/activate && pip3 install --no-cache-dir fluidfoam pandas matplotlib odfpy black'

WORKDIR /root/
RUN git clone --branch develop --recurse-submodules https://github.com/sedfoam/sedfoam.git --depth 1

WORKDIR /root
RUN /bin/bash -c 'hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -b develop && cp sedfoam/docker/Allwmakeswak swak4Foam/Libraries/Allwmake'

WORKDIR /root/swak4Foam
RUN /bin/bash -c 'shopt -s expand_aliases && ./maintainanceScripts/compileRequirements.sh'

WORKDIR /root/swak4Foam
RUN /bin/bash -c 'shopt -s expand_aliases && source /openfoam/bash.rc && export FOAM_USER_LIBBIN=$FOAM_SITE_LIBBIN && export FOAM_USER_APPBIN=$FOAM_SITE_APPBIN && export WM_NCOMPPROCS=4 && ./Allwmake'

WORKDIR /root/sedfoam
RUN /bin/bash -c 'shopt -s expand_aliases && source /openfoam/bash.rc && export FOAM_USER_LIBBIN=$FOAM_SITE_LIBBIN && ./Allwmake -prefix=group'

USER 1001

WORKDIR /home/sudofoam
RUN sudo cp -r /root/sedfoam sedfoam && sudo chown -R sudofoam: /home/sudofoam && sudo rm -rf /root/*
RUN echo "source /home/sudofoam/pyenv/bin/activate" >> .bashrc

WORKDIR /home/sudofoam
ENV HOME=/home/sudofoam
ENTRYPOINT ["/bin/bash", "-c"]
