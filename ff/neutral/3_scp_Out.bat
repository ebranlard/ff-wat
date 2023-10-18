


set simDir=/home/ebranlar/wat/ff/neutral/

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/*.outb .
scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/FF-NoWAT.out .
scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/FF-WAT.out .
