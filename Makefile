conan:
	conan install . --output-folder=build --build=missing
	cd build && cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=1

cmake_release:
	cd build && cmake --build .

cluster:
	salloc -p dev_multiple -N 2 -n 2 -t 30 # --ntasks-per-node=1 wichtig!
	module load compiler/gnu
	module load mpi/openmpi
	export OMP_NUM_THREADS=80
	mpirun --mca plm_slurm_args '--mem-per-cpu=0' build/paths # --mca plm_slurm_args '--mem-per-cpu=0' wahrscheinlich nicht wichtig

setup:
	source /etc/profile.d/modules.sh
	module load compiler/gnu
	module load mpi/openmpi
	export OMP_NUM_THREADS=80
	export MKL_NUM_THREADS=80

