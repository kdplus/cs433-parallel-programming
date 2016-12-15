#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <assert.h>

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
#endif
	return result;
}

const int TILE_DIM = 32;
const int NUM_REPS = 100;
const int SIDE = 8;
// Check errors and print GB/s
void postprocess(const float *ref, const float *res, int n, float ms)
{
	bool passed = true;
	for (int i = 0; i < n; i++)
		if (res[i] != ref[i]) {
			printf("%d %f %f\n", i, res[i], ref[i]);
			printf("%25s\n", "*** FAILED ***");
			passed = false;
			break;
		}
	if (passed)
		printf("%20.2f%20.2f\n", ms/NUM_REPS, 2 * n * sizeof(float) * 1e-6 * NUM_REPS / ms);
}

__global__ void matrixTranspose(float *_a, float *_b, const int cols, const int rows)
{
	int i = blockIdx.y * blockDim.y + threadIdx.y; // row
	int j = blockIdx.x * blockDim.x + threadIdx.x; // col
	int index_in = i*cols + j; // (i,j) from matrix A
	int index_out = j*rows + i; // transposed index
	_b[index_out] = _a[index_in];
}

__global__ void matrixTransposeShared(const float *_a, float *_b, const int cols, const int rows)
{
	__shared__ float mat[TILE_DIM][TILE_DIM];
	int bx = blockIdx.x *blockDim.x;
	int by = blockIdx.y *blockDim.y;
	int i = by + threadIdx.y; int j = bx + threadIdx.x; //input
	int ti = bx + threadIdx.y; int tj = by + threadIdx.x;

	//output
	if (i < rows && j < cols)//i < ny && j < nx
		mat[threadIdx.x][threadIdx.y] = _a[i * cols + j];
	__syncthreads(); //Wait for all data to be copied
	if (tj < cols && ti < rows)
		_b[ti * rows + tj] = mat[threadIdx.y][threadIdx.x];
}

__global__ void matrixTransposeSharedwBC(const float *_a, float *_b, const int cols, const int rows)
{
	__shared__ float mat[TILE_DIM][TILE_DIM + 1];
	int bx = blockIdx.x *blockDim.x;
	int by = blockIdx.y *blockDim.y;
	int i = by + threadIdx.y; int j = bx + threadIdx.x; //input
	int ti = bx + threadIdx.y; int tj = by + threadIdx.x;

	//output
	if (i < rows && j < cols)//i < rows && j < cols
		mat[threadIdx.x][threadIdx.y] = _a[i * cols + j];
	__syncthreads(); //Wait for all data to be copied
	if (tj < cols && ti < rows)
		_b[ti * rows + tj] = mat[threadIdx.y][threadIdx.x];
}

__global__ void matrixTransposeUnrolled(const float *_a, float *_b, const int cols, const int rows)
{
	__shared__ float mat[TILE_DIM][TILE_DIM + 1];
	int x = blockIdx.x * TILE_DIM + threadIdx.x;
	int y = blockIdx.y * TILE_DIM + threadIdx.y;
#pragma unroll
	for (int k = 0; k < TILE_DIM; k += SIDE) {
		if (x < rows && y + k < cols)
			mat[threadIdx.y + k][threadIdx.x] = _a[((y + k) * rows) + x];
	}

	__syncthreads();

	x = blockIdx.y * TILE_DIM + threadIdx.x;
	y = blockIdx.x * TILE_DIM + threadIdx.y;
#pragma unroll
	for (int k = 0; k < TILE_DIM; k += SIDE)
	{
		if (x < cols && y + k < rows)
			_b[(y + k) * cols + x] = mat[threadIdx.x][threadIdx.y + k];
	}
}

int main(int argc, char **argv){
	const int nx = 1024;
	const int ny = 1024;
	const int mem_size = nx*ny*sizeof(float);

	dim3 gridDim(nx / TILE_DIM, ny / TILE_DIM, 1);
	dim3 blockDim(TILE_DIM, TILE_DIM, 1);

	int devId = 0;
	if (argc > 1) devId = atoi(argv[1]);

	cudaDeviceProp prop;
	checkCuda(cudaGetDeviceProperties(&prop, devId));
	printf("\nDevice : %s\n", prop.name);
	printf("Matrix size: %d %d, Block size: %d %d\n",
		nx, ny, TILE_DIM, TILE_DIM);
	printf("gridDim: %d %d %d. blockDim: %d %d %d\n",
		gridDim.x, gridDim.y, gridDim.z, blockDim.x, blockDim.y, blockDim.z);

	checkCuda(cudaSetDevice(devId));

	float *h_idata = (float*)malloc(mem_size);
	float *h_tdata = (float*)malloc(mem_size);
	float *gold = (float*)malloc(mem_size);

	float *d_idata, *d_cdata, *d_tdata;
	checkCuda(cudaMalloc(&d_idata, mem_size));
	checkCuda(cudaMalloc(&d_tdata, mem_size));

	// host
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			h_idata[j*nx + i] = j*nx + i;

	// correct result for error checking
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			gold[j*nx + i] = h_idata[i*nx + j];

	// device
	checkCuda(cudaMemcpy(d_idata, h_idata, mem_size, cudaMemcpyHostToDevice));

	// events for timing
	cudaEvent_t startEvent, stopEvent;
	checkCuda(cudaEventCreate(&startEvent));
	checkCuda(cudaEventCreate(&stopEvent));
	float ms;

	// ------------
	// time kernels
	// ------------
	printf("%25s%20s%25s\n", "Method","Time(ms)", "Bandwidth (GB/s)");
	
	// ----
	// matrixTranspose 
	// ----
	printf("%25s", "matrixTranspose");
	checkCuda(cudaMemset(d_tdata, 0, mem_size));
	//matrixTranspose << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(startEvent, 0));
	for (int i = 0; i < NUM_REPS; i++)
		matrixTranspose << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(stopEvent, 0));
	checkCuda(cudaEventSynchronize(stopEvent));
	checkCuda(cudaEventElapsedTime(&ms, startEvent, stopEvent));
	checkCuda(cudaMemcpy(h_tdata, d_tdata, mem_size, cudaMemcpyDeviceToHost));
	postprocess(gold, h_tdata, nx * ny, ms);

	// ----
	// matrixTransposeShared 
	// ----
	printf("%25s", "matrixTransposeShared");
	checkCuda(cudaMemset(d_tdata, 0, mem_size));
	//matrixTransposeShared << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(startEvent, 0));
	for (int i = 0; i < NUM_REPS; i++)
		matrixTransposeShared << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(stopEvent, 0));
	checkCuda(cudaEventSynchronize(stopEvent));
	checkCuda(cudaEventElapsedTime(&ms, startEvent, stopEvent));
	checkCuda(cudaMemcpy(h_tdata, d_tdata, mem_size, cudaMemcpyDeviceToHost));
	postprocess(gold, h_tdata, nx * ny, ms);
	
	// ----
	// matrixTransposeSharedwBC 
	// ----
	printf("%25s", "matrixTransposeSharedwBC");
	checkCuda(cudaMemset(d_tdata, 0, mem_size));
	matrixTransposeSharedwBC << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(startEvent, 0));
	for (int i = 0; i < NUM_REPS; i++)
		matrixTransposeSharedwBC << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(stopEvent, 0));
	checkCuda(cudaEventSynchronize(stopEvent));
	checkCuda(cudaEventElapsedTime(&ms, startEvent, stopEvent));
	checkCuda(cudaMemcpy(h_tdata, d_tdata, mem_size, cudaMemcpyDeviceToHost));
	postprocess(gold, h_tdata, nx * ny, ms);
	
	// ----
	// matrixTransposeUnrolled 
	// ----
	dim3 blockDimUnroll(TILE_DIM, SIDE, 1);// !important
	printf("Matrix size: %d %d, Block size: %d %d\n",
	nx, ny, TILE_DIM, TILE_DIM);
	printf("gridDim: %d %d %d. blockDim: %d %d %d\n",
	gridDim.x, gridDim.y, gridDim.z, blockDim.x, blockDim.y, blockDim.z);
	printf("%25s", "matrixTransposeUnrolled");
	checkCuda(cudaMemset(d_tdata, 0, mem_size));
	//matrixTransposeUnrolled << <gridDim, blockDim >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(startEvent, 0));
	for (int i = 0; i < NUM_REPS; i++)
		matrixTransposeUnrolled << <gridDim, blockDimUnroll >> >(d_idata, d_tdata, nx, ny);
	checkCuda(cudaEventRecord(stopEvent, 0));
	checkCuda(cudaEventSynchronize(stopEvent));
	checkCuda(cudaEventElapsedTime(&ms, startEvent, stopEvent));
	checkCuda(cudaMemcpy(h_tdata, d_tdata, mem_size, cudaMemcpyDeviceToHost));
	postprocess(gold, h_tdata, nx * ny, ms);
	
error_exit:
	// cleanup
	checkCuda(cudaEventDestroy(startEvent));
	checkCuda(cudaEventDestroy(stopEvent));
	checkCuda(cudaFree(d_tdata));
	checkCuda(cudaFree(d_idata));
	free(h_idata);
	free(h_tdata);
	free(gold);
}