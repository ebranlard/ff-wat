


set simDir=/home/ebranlar/wat/ff/stable/
scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/*.outb .
:: scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/FF-NoWAT.T*.out .
:: scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/FF-WAT.T*.out .
scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/FF-NoWAT.out .
scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/FF-WAT.out .
scp -o MACs=hmac-sha2-512 ebranlar@kestrel.hpc.nrel.gov:%simDir%/FF-WAT-More.out .
