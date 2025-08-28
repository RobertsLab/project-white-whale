#!/bin/bash
# Installation script for dependencies needed by the DNA methylation data downloader

set -e

echo "Installing dependencies for DNA methylation data downloader..."

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Install SRA Toolkit
install_sra_toolkit() {
    echo "Installing SRA Toolkit..."
    
    if command_exists conda; then
        echo "Using conda to install SRA Toolkit..."
        conda install -c bioconda sra-tools
    elif command_exists apt-get; then
        echo "Using apt-get to install SRA Toolkit..."
        sudo apt-get update
        sudo apt-get install -y sra-toolkit
    elif command_exists yum; then
        echo "Using yum to install SRA Toolkit..."
        sudo yum install -y sra-toolkit
    else
        echo "Manual installation required for SRA Toolkit"
        echo "Please visit: https://github.com/ncbi/sra-tools"
        echo "Download and install according to your system"
    fi
}

# Install Entrez Direct
install_entrez_direct() {
    echo "Installing NCBI Entrez Direct..."
    
    if command_exists conda; then
        echo "Using conda to install Entrez Direct..."
        conda install -c bioconda entrez-direct
    else
        echo "Installing Entrez Direct manually..."
        cd /tmp
        curl -s https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh | bash
        echo 'export PATH=${PATH}:${HOME}/edirect' >> ~/.bashrc
        echo "Please run: source ~/.bashrc"
    fi
}

# Check what's already installed
echo "Checking existing installations..."

if command_exists fastq-dump || command_exists fasterq-dump; then
    echo "✓ SRA Toolkit is already installed"
else
    install_sra_toolkit
fi

if command_exists esearch && command_exists efetch; then
    echo "✓ Entrez Direct is already installed"
else
    install_entrez_direct
fi

# Optional tools
if command_exists parallel; then
    echo "✓ GNU parallel is available"
else
    echo "Note: GNU parallel not found. Install for faster downloads:"
    echo "  Ubuntu/Debian: sudo apt-get install parallel"
    echo "  CentOS/RHEL: sudo yum install parallel"
    echo "  macOS: brew install parallel"
fi

echo ""
echo "Installation complete!"
echo ""
echo "Test your installation with:"
echo "  fastq-dump --version"
echo "  esearch -help"
echo ""
echo "Configure SRA Toolkit with:"
echo "  vdb-config --interactive"