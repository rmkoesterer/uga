FROM ubuntu:16.04

WORKDIR /work

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
	apt-get -y --no-install-recommends install \
		locales \
		wget \
		bzip2 \
		ca-certificates && \
	sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
	locale-gen && \
	wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	/bin/bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
	rm Miniconda3-latest-Linux-x86_64.sh && \
	ln -s /opt/conda/bin/conda /usr/local/bin/conda && \
	apt-get -y remove --purge \
		wget \
		bzip2 \
		ca-certificates && \
	apt-get -y autoremove && \
	rm -rf /var/lib/apt/lists/*

# Set the locale
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8  
ENV LANGUAGE=en_US:en  
ENV LC_ALL=en_US.UTF-8     

RUN conda config --set auto_update_conda false && \
	echo "name: uga" > /work/environment.yml && \
	echo "channels:" >> /work/environment.yml && \
	echo "- conda-forge" >> /work/environment.yml && \
	echo "- bioconda" >> /work/environment.yml && \
	echo "- r" >> /work/environment.yml && \
	echo "- defaults" >> /work/environment.yml && \
	echo "dependencies:" >> /work/environment.yml && \
	echo "- python=3.7.3" >> /work/environment.yml && \
	echo "- numpy=1.16.4" >> /work/environment.yml && \
	echo "- nomkl=3.0" >> /work/environment.yml && \
	echo "- pandas=0.23.4" >> /work/environment.yml && \
	echo "- psutil=5.6.3" >> /work/environment.yml && \
	echo "- pysam=0.15.2" >> /work/environment.yml && \
	echo "- scipy=1.3.0" >> /work/environment.yml && \
	echo "- setuptools=41.0.1" >> /work/environment.yml && \
	echo "- cython=0.29.12" >> /work/environment.yml && \
	echo "- biopython=1.74" >> /work/environment.yml && \
	echo "- progressbar=2.5" >> /work/environment.yml && \
	echo "- r-base=3.5.1" >> /work/environment.yml && \
	echo "- rpy2=2.9.4" >> /work/environment.yml && \
	conda env create -f /work/environment.yml && \
	conda clean -afy

RUN echo "conda activate uga" >> ~/.bashrc
ENV PATH /opt/conda/envs/uga/bin:/opt/conda/bin:$PATH

RUN R -e 'install.packages(c("kinship2", "geepack", "lme4", "lmertest", "pbkrtest", "seqMeta", "RColorBrewer", "ggplot2","R.utils"), repos="http://cran.us.r-project.org", dependencies=TRUE)'

RUN apt-get update && \
	apt-get -y --no-install-recommends install \
		gcc \
		python-dev && \
	pip install uga && \
	ln -s /opt/conda/envs/uga/bin/python /usr/local/bin/python && \
	chmod 775 /opt/conda/envs/uga/lib/python2.7/site-packages/uga/Qsub.py && \
	chmod -R 777 /work && \
	apt-get -y remove --purge \
		gcc \
		python-dev && \
	apt-get -y autoremove && \
	rm -rf /var/lib/apt/lists/*

RUN echo "#!/bin/bash" > /work/.uga_init && \
	echo ". activate uga" >> /work/.uga_init && \
	chmod 775 /work/.uga_init

ENTRYPOINT ["/work/.uga_init"]
