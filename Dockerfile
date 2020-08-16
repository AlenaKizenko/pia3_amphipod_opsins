FROM continuumio/miniconda3

MAINTAINER Maria Firulyova <mmfiruleva@gmail.com>

ADD smk.yml /tmp/smk.yml

RUN apt-get update && \
    apt-get install -y libfftw3-double3 && \
    conda env update -n base --file /tmp/smk.yml && \
    mkdir -p /home/PIA3

ADD PIA3 /home/PIA3

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

WORKDIR /home/PIA3
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]