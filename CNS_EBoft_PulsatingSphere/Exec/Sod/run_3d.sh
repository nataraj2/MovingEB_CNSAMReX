rm -rf CNS_EBoft3d.gnu.MPI.ex
rm -rf plt*
make
/usr/local/Cellar/mpich/3.3.1/bin/mpiexec -np 8 ./CNS_EBoft3d.gnu.MPI.ex inputs

