FROM openfoamplus/of_v2006_centos73
LABEL maintainer "cyrille.bonamy@univ-grenoble-alpes.fr"
ARG WM_NCOMPPROCS=10

RUN yum update -y \
  && yum install -y \
  git python3 python3-pip unzip \
  mercurial readline-devel \
  texlive dvipng python3-tkinter \
  && yum clean all
RUN pip3 install fluidfoam

USER 1000

WORKDIR /home/ofuser
RUN git clone --recurse-submodules https://github.com/sedfoam/sedfoam.git

WORKDIR /home/ofuser
RUN /bin/bash -c "hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam -r v2020.06 && cp sedfoam/docker/Allwmakeswak swak4Foam/Libraries/Allwmake"

WORKDIR /home/ofuser/swak4Foam
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./maintainanceScripts/compileRequirements.sh"
WORKDIR /home/ofuser/swak4Foam
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./Allwmake"
#RUN /bin/bash -c "source /opt/openfoam/5.0/OpenFOAM-5.0/etc/bashrc && ./maintainanceScripts/copySwakFilesToSite.sh"

WORKDIR /home/ofuser
RUN rm -rf swak4Foam

WORKDIR /home/ofuser/sedfoam
RUN /bin/bash -c "shopt -s expand_aliases && source /opt/OpenFOAM/setImage_v2006.sh && ./Allwmake && ./Allwclean"
#RUN find /home/ofuser/sedfoam -name "*.o" -exec rm -f {} \;

WORKDIR /home/ofuser
ENV HOME /home/ofuser
ENTRYPOINT ["/bin/bash", "-c"]
