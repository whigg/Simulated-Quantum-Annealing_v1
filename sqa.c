#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
MC_stpes:
M: number of the trotters
N: number of the spin
--------------------------------
SQA algo_example
1 for(int t = 0; t<T; t++) {
2 	for(int m=0; m<M; m++) {
3 		for(i=0; i<N; i++) {
4 			local_field[i] = 0;
5 			for(j=0; j<N; j++) {
6 				//Compute energy
7 				local_field[i] += spin(m,j) * J(i,j);
8 			}
9 			//Compute energy considering the two neighboring Trotter slices;
10 			local_field[i] -= spin(m-1,i) * K;
11 			local_field[i] += spin(m+1,i) * K;
12 			if( exp(-1.0*local_field[i]/T) > rand_num )
13 				spin[i] = !spin[i];
14 			}
15 		}
16   }
17 }

*/
#define MC_steps 100
#define M 64
#define N 2048 



void usage(){
	printf("Usage:\n");
	printf("       ./sqa [spin configuration]\n");
	exit(0);
}

int main (int argc, char *argv[]){
	if (argc != 2)
		usage();

	//x,y <- coordinate, and weight
	//beta = -1/(k_b*T), T is temperature, here means MC_steps, use metrapolesis
	double k_b = 1.381e-23;
	int x, y, idx;
	int wei = N*N;
	int results = 0;
	int beta;
	
	// couplings, spin, single Hamiltonion energy
	char *J;
	char *spin;
	int *local_field;
	J = (char*)malloc(N*N*sizeof(char));
	spin = (char*)malloc(M*N*sizeof(char));
	local_field = (int*)malloc(M*N*sizeof(int));

	//initialize: couplings, sgl_Jam zeroing, spin random num, Hamil_init
	for(int i = 0; i < wei; i++){
		J[i] = 0;
	}
	// giving spin random values
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            spin[m * N + n] = (rand() / (float) RAND_MAX) > .5 ? 1 : -1;
        }
    }
    
	//read file
	FILE *file = fopen(argv[1], "r");
	if(file == NULL){
		printf("error!\n");
		return 1;
	}
	fscanf(file, "%d", &wei);
	while(!feof(file)){
		fscanf(file, "%d%d%d",  &x, &y, &wei);
		if (x != y ){
			J [x * N + y] = wei;
			J [y * N + x] = wei;
		}
	}
	fclose(file);

	//starts, only run one time
	for (int t = 0; t < MC_steps; t++){
		for (int m = 0; m < M; m++){
			for (int i = 0; i < N; i++){
				//one spin flip computation
				//pseudo: local_field[i] = 0;
				idx = m*N+i;
				local_field[idx] = 0;
				//compute local field 
				for(int j = 0; j < N; j++){
					//Compute energy
					//pseudo: local_field[i] += spin(m,j) * J(i,j);
					local_field[idx] += spin[m*N+j] * J[i*N+j];
				}
				//Compute energy considering the two neighboring Trotter slices == compute transverse-field
				//local_field += Jtrans * (spin(m+1,i)-spin(m-1,i))
				//energy_diff = local_field*...
				////---------------這裡要補：如果在最上層或者是最下層的時候要怎辦
				local_field[idx] -= spin[(m-1)*N+i]*k_b;
				local_field[idx] += spin[(m+1)*N+i]*k_b;
				//pseudo: local_field[i] -= spin(m-1,i) * K;
			 	//pseudo: local_field[i] += spin(m+1,i) * K;
	
				//spin update 
				if( exp(-1.0*local_field[idx]/MC_steps) > (rand()/(float)RAND_MAX) ){
 					spin[idx] = !spin[idx];
				}
				//if(exp(-energy_diff/....)>rand_num){
				//	spin[i] = negate spin[i]
				//}
			}
		}
	}	
	int E = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            E += -spin[i]*spin[j]*J[i*N+j];
    results = E;


	FILE *output;
	output = fopen("sqa_ful_cc.txt", "w");
	fprintf(output, "%d\n", results);
	fclose(output);
	free(J);
	free(spin);
	free(local_field);

	return 0;
}
