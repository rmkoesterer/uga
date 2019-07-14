FROM python:2.7.16-slim-stretch

# all main work in work
WORKDIR /work

# open /work directory permissions
RUN chmod -R 777 /work

# https support, libcurl, locales
RUN apt-get update && \
	apt-get -y install apt-utils apt-transport-https libcurl4-openssl-dev locales gnupg2 software-properties-common

# Set the locale
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8  
ENV LANGUAGE=en_US:en  
ENV LC_ALL=en_US.UTF-8     

# install debian package for R 3.4
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' && \
	apt-get update && \
	apt-get -y install r-base && \
	ln -s /usr/bin/R /usr/local/bin/R && \
	ln -s /usr/bin/Rscript /usr/local/bin/Rscript

# install R packages
RUN R -e 'install.packages(c("kinship2", "geepack", "lme4", "lmertest", "pbkrtest", "seqMeta", "RColorBrewer", "ggplot2","R.utils"), repos="http://cran.us.r-project.org", dependencies=TRUE)'

RUN pip install demandimport && \
	pip install cython && \
	pip install numpy && \
	pip install pandas && \
	pip install psutil && \
	pip install pysam && \
	pip install rpy2==2.8.6 && \
	pip install scipy && \
	pip install setuptools && \
	pip install biopython && \
	pip install progressbar

# install uga
#RUN pip install uga
