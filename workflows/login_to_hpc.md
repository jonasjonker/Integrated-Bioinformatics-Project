# How to: login to the HPC

### Generate an id_rsa key-pair
If you have never loged in to the HPC you need to upload a ssh key. The HPC
requires bigger keys than the default, so you might need to generate a new one
using the following command:

```bash
ssh-keygen -t rsa -b 4096
```
This program will ask how you wnat to name your key. If you already have rsa
keys on your machine you might need to give a non default name or different path.
for example: `~/.ssh/id_rsa_4096`.

### upload the public key
Login on [vscentrum](https://account.vscentrum.be/) with your KU Leuven credentials.

Go to the tab _Edit Acount_, and _upload your **public** key_

It might take 15 minutes before the HPC recognizes your new key.

### login
```bash
ssh <your_vsc_id>@login.hpc.kuleuven.be
```
