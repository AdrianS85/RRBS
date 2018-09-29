BootStrap: docker
From: ubuntu:18.04
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

#BUILD: sudo singularity build --writable RRBS_Singularity_New.simg ./Singularity2
#SANDBOX: sudo singularity build -s sandbox1.simg ./Singularity2
#ENTER: sudo singularity shell --cleanenv --bind ./Bind:/ --workdir /Analysis --writable RRBS_Singularity_New.simg

#File folder needs to be structured as such: 
#Create directory which will be binded with singularity on the root level (for example Bind: --bind ./Bind:/). This will take care of temporary files being written in singularity-based folders that would not have access to disk space.
#Create directory inside this directory, which will hold all the data (for example Analysis)
#Hence in a single directory we will have .simg file and Bind/Analysis directory. The Bind/Analysis directory must contain:
#1) Files to be analyzed (.fastq format, .gz packed)
#2) Folder "Genome" with Bismark index "Bisulfite_Genome" and reference genome in .fa format
#3) nextflow workflow file

%environment
        export PATH=/samtools-1.9:$PATH
	export PATH=/Bismark-master:$PATH
	export PATH=/FastQC:$PATH
	export PATH=/TrimGalore-master:$PATH
	export PATH=/BamQC-master/bin:$PATH
	export PATH=/nudup-master:$PATH
	export LANGUAGE=en_US.UTF-8
        export LANG=en_US.UTF-8
        export LC_ALL=en_US.UTF-8

%post
        apt update
        apt -y install vim wget perl unzip default-jdk bowtie2 python-pip libcurl3 libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev ant language-pack-en
	
        cd /
        
        wget https://github.com/nugentechnologies/NuMetRRBS/blob/master/strip_bismark_sam.sh
        mv strip_bismark_sam.sh /usr/bin
        
        wget https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py
        
        
        wget -qO- https://get.nextflow.io | bash ####v0.32.0.4897
        mv nextflow /usr/bin
	/usr/bin/nextflow -v ####nextflow downloads additional dependencies here

        wget https://github.com/FelixKrueger/Bismark/archive/master.zip ####v0.20.0
        unzip master.zip
        rm master.zip
        
        wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
        unzip fastqc_v0.11.7.zip
        chmod 755 FastQC/fastqc
        rm fastqc_v0.11.7.zip
        
        wget https://github.com/FelixKrueger/TrimGalore/archive/master.zip ####v0.5.0_dev
        unzip master.zip
        rm master.zip
        
        wget https://github.com/s-andrews/BamQC/archive/master.zip ####v0.1.25_devel
        unzip master.zip
	cd /BamQC-master
	ant
	chmod 755 bin/bamqc
	cd ..
        rm master.zip
        
        wget https://github.com/nugentechnologies/nudup/archive/master.zip
        unzip master.zip
        rm master.zip
        
        wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
        tar -xjvf samtools-1.9.tar.bz2
        cd samtools-1.9
        ./configure
        make
        make install
        cd ..
        rm samtools-1.9.tar.bz2
        
	pip install multiqc cutadapt

	#### multiqc.v1.6, cutadapt.v1.18, bowtie2.v2.3.4.1.64bit, openjdk.v10.0.2, Python 2.7.15rc1
