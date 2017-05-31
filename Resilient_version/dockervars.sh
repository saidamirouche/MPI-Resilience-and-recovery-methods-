#/bin/sh
#
# Load some alias to compile and run with the docker mpicc/mpif90/mpirun
#

case $1 in
    load)
        alias make='docker run -v $PWD:$PWD -w $PWD ftmpi/ulfm:1.1 make'
        alias mpirun='docker run -v $(pwd):$(pwd) -w $(pwd) ftmpi/ulfm:1.1 mpirun'
        ;;
    unload)
        unalias make
        unalias mpirun
        ;;
    *)
        echo ". $0 load: alias make and mpirun in the current shell"
        echo ". $0 unload: remove aliases from the current shell"
esac


