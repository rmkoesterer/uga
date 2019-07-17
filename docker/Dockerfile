FROM debian:stretch-slim

# all main work in work
WORKDIR /work

ENV DEBIAN_FRONTEND noninteractive

# https support, libcurl, locales
RUN apt-get update && \
	apt-get -y --no-install-recommends install \
		apt-utils \
		apt-transport-https \
		libcurl4-openssl-dev \
		locales \
		gnupg2 \
		software-properties-common \
		dirmngr \
		wget \
		bzip2 \
		ca-certificates \
		libglib2.0-0 \
		libxext6 \
		libsm6 \
		libxrender1 && \
	sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && locale-gen && \
	echo "deb http://cran.rstudio.com/bin/linux/debian stretch-cran34/" >> /etc/apt/sources.list && \
	apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	apt-get update && \
	apt-get -y install r-base && \
	ln -s /usr/bin/R /usr/local/bin/R && \
	ln -s /usr/bin/Rscript /usr/local/bin/Rscript && \
	R -e 'install.packages(c("kinship2", "geepack", "lme4", "lmertest", "pbkrtest", "seqMeta", "RColorBrewer", "ggplot2","R.utils"), repos="http://cran.us.r-project.org", dependencies=TRUE)' && \
	wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
	/bin/bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/conda && \
	rm Miniconda2-latest-Linux-x86_64.sh && \
	ln -s /opt/conda/bin/conda /usr/local/bin/conda && \
	apt-get -y remove --purge \
		apt-transport-https \
		libcurl4-openssl-dev \
		software-properties-common \
		wget \
		ca-certificates \
		libxext6 \
		libsm6 \
		libxrender1 \
		apt-utils \
		dirmngr \
		gnupg2 && \
	apt-get -y autoremove && \
	rm -rf /var/lib/apt/lists/*

# update locale
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# find all versions that are compatible with python2.7
RUN conda config --set auto_update_conda false && \
	echo "name: uga" > /work/environment.yml && \
	echo "channels:" >> /work/environment.yml && \
	echo "- r" >> /work/environment.yml && \
	echo "- defaults" >> /work/environment.yml && \
	echo "- bioconda" >> /work/environment.yml && \
	echo "- mittner" >> /work/environment.yml && \
	echo "- https://conda.anaconda.org/rmkoesterer" >> /work/environment.yml && \
	echo "dependencies:" >> /work/environment.yml && \
	echo "- python=2.7" >> /work/environment.yml && \
	echo "- numpy" >> /work/environment.yml && \
	echo "- nomkl" >> /work/environment.yml && \
	echo "- pandas" >> /work/environment.yml && \
	echo "- psutil" >> /work/environment.yml && \
	echo "- pysam=0.15.2" >> /work/environment.yml && \
	echo "- rpy2" >> /work/environment.yml && \
	echo "- scipy" >> /work/environment.yml && \
	echo "- setuptools" >> /work/environment.yml && \
	echo "- cython" >> /work/environment.yml && \
	echo "- biopython" >> /work/environment.yml && \
	echo "- progressbar" >> /work/environment.yml && \
	conda env create -f /work/environment.yml && \
	conda clean -afy

ENV PATH /opt/conda/envs/uga/bin:/opt/conda/bin:$PATH

RUN echo "source activate uga" > ~/.bashrc && \
	pip install uga==2.0.4 && \
	ln -s /opt/conda/envs/uga/bin/python /usr/local/bin/python && \
	chmod 775 /opt/conda/envs/uga/lib/python2.7/site-packages/uga/Qsub.py && \
	chmod -R 777 /work
