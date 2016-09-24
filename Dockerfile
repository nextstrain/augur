# nextstrain-augur dockerfile

FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy
RUN apt-get install -y libpng-dev libfreetype6-dev pkg-config
RUN apt-get install -y libatlas-base-dev

# mafft
RUN apt-get install -y mafft

# fasttree
RUN apt-get install -y fasttree

# raxml
RUN apt-get install -y raxml
RUN cp /usr/bin/raxmlHPC /usr/bin/raxml


# TreeTime
RUN git clone https://github.com/neherlab/treetime.git /TreeTime
WORKDIR /TreeTime
RUN python setup.py install

# augur
RUN git clone https://github.com/nextstrain/augur.git /augur
WORKDIR /augur

# python modules
RUN pip install biopython==1.64
RUN pip install cvxopt --user
RUN pip install DendroPy==3.12.0
RUN pip install boto==2.38.0
RUN pip install matplotlib==1.5.1
RUN pip install pandas==0.16.2
RUN pip install seaborn==0.6.0

# This should be used once the requirements.txt file is added to master branch
# so that only one file needs to be modified to change dependencies.
# RUN pip install -r requirements.txt

# update
ADD http://www.timeapi.org/utc/now /tmp/bustcache
RUN git pull

# default process
CMD /bin/bash