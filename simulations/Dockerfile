FROM condaforge/mambaforge:23.3.1-1
LABEL authors="Kaido Lepik" \
      description="Docker image containing all requirements for the mediation project simulations"

# Install R dependencies with conda/mamba
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -a
# Put the installed dependencies to PATH
ENV PATH /opt/conda/envs/environment/bin:$PATH

# Install plink1.9 (2021-06-06) to /usr/local/bin
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip \
 && unzip plink_linux_x86_64_20210606.zip -d Plink1.9 \
 && mv Plink1.9/plink /usr/local/bin/plink1.9 \
 && ln -s /usr/local/bin/plink1.9 /usr/local/bin/plink \
 && rm -r Plink1.9 plink_linux_x86_64_20210606.zip

# Install TwoSampleMR
RUN R -e "remotes::install_github('qingyuanzhao/mr.raps'); if (!library(mr.raps, logical.return = TRUE)) quit(status = 10)" \
 && R -e "remotes::install_github('MRCIEU/TwoSampleMR@v0.5.7'); if (!library(TwoSampleMR, logical.return = TRUE)) quit(status = 10)"
