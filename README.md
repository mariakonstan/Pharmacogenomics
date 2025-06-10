# Pharmacogenomics
WGS scanning tool

## Requirements
- Linux OS or Windows WSL
- Python 3.10+
- samtools 1.22

### Installing samtools
- Install Updates and Required Packages
```bash
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install make
sudo apt-get install bzip2
sudo apt-get install libbz2-dev
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev 
sudo apt-get install libncursesw5-dev
sudo apt-get install liblzma-dev
```

- Install HTSLIB
```bash
cd /usr/bin
sudo wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2
sudo tar -vxjf htslib-1.22.tar.bz2
cd htslib-1.22
sudo make
```

- Install SAMTOOLS
```bash
cd /usr/bin
sudo wget https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2
sudo tar -vxjf samtools-1.22.tar.bz2
cd samtools-1.22
sudo make
```

- Install BCFTools
```bash
cd /usr/bin
sudo wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
sudo tar -vxjf bcftools-1.22.tar.bz2
cd bcftools-1.22
sudo make
```

- Export the binaries to `PATH`
```bash
vi ~/.profile
```

Paste these lines at the bottom of the file...
```bash
# Add bioinformatics tools to PATH
export PATH="$PATH:/usr/bin/bcftools-1.22:/usr/bin/samtools-1.22:/usr/bin/htslib-1.22"
```

Save and close the file by typing `wq` and Enter.

Apply the changes to the current session.
```bash
source ~/.profile
```

## Setup
- Create a virtual environment
```bash
python3 -m venv venv
```

- Activate the virtual environment
```bash
source venv/bin/activate
```

- Install packages
```bash
pip install -r requirements.txt
```

## Run
- Make any adjustments to the `project/settings.py` (if necessary).
- Execute...
```bash
python3 scan.py
```

## Post-Processing
- Generate csv files and plots:
```bash
python3 process.py
```