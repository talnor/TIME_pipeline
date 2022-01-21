FROM python:2

RUN apt-get update && \
  apt-get install -yq make git  && \
  apt-get clean

RUN pip install --upgrade pip
RUN pip install 'numpy==1.16.6' && \
pip install 'biopython==1.68' 'pyfastaq==3.17.0'

RUN curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O && \
    bash Miniconda3-4.7.12.1-Linux-x86_64.sh -bf -p /usr/miniconda3/ && \
    rm Miniconda3-4.7.12.1-Linux-x86_64.sh && \
    /usr/miniconda3/bin/conda clean -tipsy && \
    ln -s /usr/miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /usr/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
COPY requirements.yml /requirements.yml
RUN /usr/miniconda3/bin/conda env create -f /requirements.yml && \
/usr/miniconda3/bin/conda clean -a

ENV PATH=/usr/miniconda3/envs/time_analysis/bin/:$PATH

CMD [ "/bin/bash" ]
