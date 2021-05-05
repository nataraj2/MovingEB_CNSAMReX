rm -rf CNS3d.gnu.MPI.ex
rm -rf plt*
make
/usr/local/Cellar/mpich/3.3.1/bin/mpiexec -np 8 ./CNS3d.gnu.MPI.ex inputs

