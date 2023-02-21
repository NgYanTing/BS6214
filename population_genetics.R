#### Code records for assignment 3 ####
# sratoolkit.3.0.0-mac64/bin already in ~/.bashrc
source ~/.bashrc

prefetch SRR6156276 
prefetch SRR6156277 --max-size 400000000000

# To login to NSCC
ssh -o HostKeyAlgorithms=+ssh-rsa -o PubkeyAcceptedKeyTypes=+ssh-rsa ngya0020@ntu.nscc.sg
