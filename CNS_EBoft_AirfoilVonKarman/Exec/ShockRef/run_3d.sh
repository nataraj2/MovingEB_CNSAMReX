rm -rf CNS_EBoft_actual3d.gnu.MPI.ex
rm -rf chk*
rm -rf plt*
make
/usr/local/Cellar/mpich/3.3.1/bin/mpiexec -np 16 ./CNS_EBoft_actual3d.gnu.MPI.ex inputs

