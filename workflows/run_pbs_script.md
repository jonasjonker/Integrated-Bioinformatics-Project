# How to: run .pbs scripts on the HPC

### Login
```bash
~$ ssh <your_vsc_id>@login.hpc.kuleuven.be
```

### Setup your HPC account for _singularity_ use
`$HOME` is to small for singularity cache. Use `$VSC_SCRATCH` instead by adding
the following lines to your `.bashrc` _dotfile_.
```bash
if [[ ! -d "$SINGULARITY_TMPDIR" ]]; then
    mkdir -p $SINGULARITY_TMPDIR
fi

if [[ ! -d "$SINGULARITY_CACHEDIR" ]]; then
    mkdir -p $SINGULARITY_CACHEDIR
fi
export SINGULARITY_TMPDIR=$VSC_SCRATCH/singularity_tmp
export SINGULARITY_CACHEDIR=$VSC_SCRATCH/singularity_cache
```

### Create an _singularity_ image
There are multiple ways to create _singularity_ images. The easiest is to
use a docker images as base.

First, make a directory in $VSC_DATA to store you images. And make a symlink to
it.
```bash
~$ mkdir ${VSC_DATA}/images
~$ ln -s $(VSC_DATA)/images
```

Then build an image. (The example image below is very small. Just try it)
```bash
#  singularity build <name> docker://<docker_image>
~$ singularity build images/alpine.sif docker://alpine:latest
~$ ls images
```

### Use a _singularity_ image
Execute a command (like a script) within a singularity container.
```bash
~$ singularity exec images/alpine.sif cat /etc/os-release
NAME="Alpine Linux"
ID=alpine
VERSION_ID=3.12.0
PRETTY_NAME="Alpine Linux v3.12"
HOME_URL="https://alpinelinux.org/"
BUG_REPORT_URL="https://bugs.alpinelinux.org/"
```

you can also open a terminal within a container.
```bash
~$ singularity shell alpine.sif
singularity> ...
```

Or execute the Images predefined commands (if any)
```bash
~$ singularity run alpine.sif
# for this specific image, it just opens a terminal.
singularity> ...
```

### Share a _singularity_ image
You can share singularity images (or any file in the hpc) by changing the file
 permissions
```bash
# gives everybody read and execute permission, but only you write permission.
~$ chmod 755 images/alpine.sif
```
