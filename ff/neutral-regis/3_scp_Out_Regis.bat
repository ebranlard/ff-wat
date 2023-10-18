


set simDir=/home/ebranlar/wat/ff/neutral-regis/regis/Seed_0/

scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/*.out .
scp -o MACs=hmac-sha2-512 ebranlar@eagle.hpc.nrel.gov:%simDir%/*.dbg .
