#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <windows.h>
#include <direct.h>


//constants
#define NUM_INTERATIONS 200									//Number of iterations for the reevaluation of our model parameters
#define M 32 												//The number of distinct observational symbols per state
#define N 5 												//The number of states in the model
#define P 12												//Number of features in the cepstral coefficients
#define LIMIT 5000											//The max absolute value of a sample
#define CODEBOOKSIZE 32										//size of the codebook
#define PI 3.142857142857									//22/7
#define FRAME_SIZE 320										//size of a frame
#define FRAMESMAXSIZE 100									//max num of frames to consider in a speech
#define SAMPLESMAXSIZE 32000								//max num of samples to consider in a speech
const int MAX_T = 150; 										//The maximum number of observations in the sequence

//variables
int T;
long double dcshift, multfact, absmax;
long double const zerosubstitute = 1e-30;
long int maxsamplesize = 0, samplesize = 0, numframes = 0;
long double maximumprob = 0;
int answer = 0;
int O[MAX_T+1];	
int Q[MAX_T+1];	
long double probgivenOandLambda = 0, XprobgivenOandLambda = -1;
long double Alpha[MAX_T+1][N+1];
long double Beta[MAX_T+1][N+1];
long double Gamma[MAX_T+1][N+1];
long double Delta[MAX_T+1][N+1];
int Shi[MAX_T+1][N+1]; 
long double Xi[MAX_T+1][N+1][N+1];
long double codeBook[CODEBOOKSIZE][P];
long int sample[1000000];
long double inputsamples[SAMPLESMAXSIZE];
long double Ai[P+1], Ri[P+1], Ci[P+1];
long double A[N+1][N+1] = {0};
long double B[N+1][M+1] = {0};
long double Pi[N+1] = {0};
long double Anew[N+1][N+1] = {0};
long double Bnew[N+1][M+1] = {0};
long double Pinew[N+1] = {0};
long double Aaverage[N+1][N+1] = {0};
long double Baverage[N+1][M+1] = {0};
long double Piaverage[N+1] = {0};
char* locationA = "A.txt";
char* locationB = "B.txt";
char* locationPI = "PI.txt";
int count = 1, train = 0;
long double probOLambda = 0;


//weights for measuring tokhura distances
double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//forward procedure to calculate the prob of observing a observation sequence given a model 
void forwardProcedure(){
	int i , j , t;
	long double sum ;
	int index = O[1];
	probOLambda = 0;

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}
	
	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	for(i=1;i<=N;i++){
		probOLambda = probOLambda + Alpha[T][i];
	}
}

//It finds the digit of the model that has the highest value of the prob P(O|(A,B,Pi)) 
void forwardProcedurePredict(int iteration){
	int i , j , t;
	long double sum ;
	int index = O[1];
	probOLambda = 0;

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}
	
	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	for(i=1;i<=N;i++){
		probOLambda = probOLambda + Alpha[T][i];
	}
	
	if(probOLambda > maximumprob){
		maximumprob = probOLambda;
	 	answer = iteration;
	}
}

//backward procedure to calculate beta
void backwardProcedure(){
	int i , j , t;
	long double sum;
	int index = 0;
	for(i=1;i<=N;i++){
		Beta[T][i] = 1.0;
	}
	for(t=T-1;t>=1;t--){
		index = O[t+1];
		for(i=1;i<=N;i++){
			sum = 0;
			for(j=1;j<=N;j++){
				sum = sum + B[j][index]*A[i][j]*Beta[t+1][j];
			}
			Beta[t][i]=sum;
		}
	}
}

//procedure to find Gamma
void findGamma(){
	for(int t=1;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double summation=0;
			for(int k=1;k<=N;k++){
				summation += Alpha[t][k] * Beta[t][k];
			}
			Gamma[t][j]=(Alpha[t][j] * Beta[t][j])/summation;
		}
	}
}

//function to update the parameters of our model with newly evaluated parameters
void loadModel(){
	int i, j;
	for(i=1;i<=N;i++){
		Pi[i]=Pinew[i];
	}
	
	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			A[i][j]= Anew[i][j];
		}
	}
	for(i=1;i<=N;i++){
		for(j=1;j<=M;j++){
			B[i][j] = Bnew[i][j];
		}
	}
}

//the function assigns a very small number in place of zeroes in the elements of B
void applyZeroSubstituteToB(){
	int i, j;
	long double diff;
	long double max;
	int max_i=0;
	int max_j=0;
	for (i = 1; i <= N; i++){
		diff = 0;
		max = 0;
		for (j = 1; j <= M; j++){
			if (Bnew[i][j] > max){
				max = Bnew[i][j];
				max_i = i;
				max_j = j;
			}
			if (Bnew[i][j] < zerosubstitute){
				diff += Bnew[i][j] - zerosubstitute;
				Bnew[i][j] = zerosubstitute;
			}
		}
		Bnew[max_i][max_j] = max;
	}
}

//Reevaluate the parameters of the model
void recalculateModel(){
	int i, j, k, t;
	long double sum1=0 , sum2 =0;
	for(i=1;i<=N;i++){
		Pinew[i] = Gamma[1][i];
	}
	
	for(i = 1; i<=N; i++){
		for(j = 1; j <= N; j++){
			long double t1 = 0, t2 = 0;
			for(t = 1; t <= T-1; t++){
				t1 += Xi[t][i][j];
				t2 += Gamma[t][i];
			}
			Anew[i][j] = t1/t2;
		}
	}
	
	for(j=1;j<=N;j++){
		int count=0;
		long double max=0;
		int ind_j=0, ind_k=0;
		
		for(k=1;k<=M;k++){
			sum1 =0 , sum2 =0;
			for(t=1;t<T;t++){
				sum1 = sum1 + Gamma[t][j];
				if(O[t]==k){
					sum2 = sum2 + Gamma[t][j];				
				}
			}
			Bnew[j][k] = sum2/sum1;
			
			
			if(Bnew[j][k]>max){
				max=Bnew[j][k];
				ind_j = j;
				ind_k = k;
			}
			
			
			if(Bnew[j][k] == 0){
				Bnew[j][k]=zerosubstitute;
				count++;
			}
		}
		Bnew[ind_j][ind_k] = max - count*zerosubstitute;
	}
	
	loadModel();
}

//Finds the matrix Xi
void findXi(){

	int i , j , t;
	long double summation[FRAMESMAXSIZE + 1];

	for(t=1;t<=T;t++){
		
		summation[t] = 0;
		for(i=1;i<=N;i++){
			for(j=1;j<=N;j++){
				summation[t] += Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
			}
		}

		for(i=1;i<=N;i++){
			long double x;
			for(j=1;j<=N;j++){
				x = Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
				Xi[t][i][j]= x/summation[t];
			}
		}
	}
}

//Viterbi algorithm for solving HMM Problem 2
void viterbiAlgorithm(){
	
    for(int i=1; i<=N; i++){
        Delta[1][i] = Pi[i] * B[i][O[1]];
        Shi[1][i] = 0;
    }

	
	for(int j=1; j<=N; j++){
		for(int t=2; t<=T; t++){
            long double max = 0, ti = 0;
            int ind = 0;
            
            for(int i=1; i<=N; i++){
                ti = Delta[t-1][i] * A[i][j];
                if(ti > max){
					max = ti;
					ind = i;
				}
            }

            Delta[t][j] = max * B[j][O[t]];
			Shi[t][j] = ind;
        }
    }


    long double max = 0;
    for(int i=1; i<=N; i++){
        if(Delta[T][i] > max) {
			max = Delta[T][i];
			Q[T] = i;
		}

        probgivenOandLambda = max;
    }


    for(int t = T-1; t>0; t--){
        Q[t] = Shi[t+1][Q[t+1]];
    }
}

//Read matrix A from file
bool readA(char *filename){
    FILE* fp = fopen(filename, "r");												
	if(!fp)
	{
		printf("Couldn't read file!\n");
		return false;
	}
	while(!feof(fp)){
		for(int i = 1; i <= N; i++){
			for(int j = 1; j <= N; j++)
				fscanf(fp, "%Lf", &A[i][j]);
		}
	}

	fclose(fp);
	return true;
}

//Read matrix B from file
bool readB(char *filename){
	FILE* fp = fopen(filename, "r");												

	while(!feof(fp)){
		for(int i = 1; i <= N; i++){
			for(int j = 1; j <= M; j++)
				fscanf(fp, "%Lf", &B[i][j]);
		}
	}

	fclose(fp);
	return true;
}

//Read Pi from file
bool readPi(char* filename){
	FILE *fp = fopen(filename, "r");

	while(!feof(fp)){
		for(int i = 1; i <= N; i++)
			fscanf(fp, "%Lf", &Pi[i]);
	}

	fclose(fp);
	return true;
}

//Function to erase models' parameters
void deleteModels(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
			Aaverage[i][j] = 0;
			Anew[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
			Baverage[i][j] = 0;
			Bnew[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
		Pinew[i] = 0;
		Piaverage[i] = 0;
	}
}

//fnction to erase parameters of the global model
void deleteTheModel(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
	}
}

//function to erase parameters of the avg global model
void deleteAvgModel(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			Aaverage[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			Baverage[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Piaverage[i] = 0;
	}
}

//function to read saved model from disk
void readAvgModel(int digit, char* user){
	
	char filename[100];
	sprintf(filename, "output/savedmodels/%s/digit_%d_A.txt", user, digit);
	readA(filename);

	sprintf(filename, "output/savedmodels/%s/digit_%d_B.txt", user, digit);
	readB(filename);

	sprintf(filename, "output/savedmodels/%s/digit_%d_PI.txt", user, digit);
	readPi(filename);
}

//initialize the model
void initializeModel(int digit, int seq, char* user, char *filename = "skeletalModel"){

	if(filename == "skeletalModel"){
		readA(locationA);
		readB(locationB);
		readPi(locationPI);
	}else if(filename  == "avg"){
		readAvgModel(digit, user);
	}else if(filename == "init"){
		readA(locationA);
		readB(locationB);
		readPi(locationPI);
	}
}

//adds parameters of a model to the average model
void addToAvgModel(){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			Aaverage[i][j] += A[i][j];
		}
	}
	for (i = 1; i <= N; i++){
		Piaverage[i] += Pi[i];
	}
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= M; j++){
			Baverage[i][j] += B[i][j];
		}
	}
}

//Update the matrix A in file
void updateA(char filename[]){
	FILE *fp = fopen(filename, "w");
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			
			fprintf(fp, "%Le   ", A[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//update the matrix B in file
void updateB(char filename[]){
	FILE *fp = fopen(filename, "w");
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= M; j++){
			
			fprintf(fp, "%Le   ", B[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//update Pi in file
void updatePI(char filename[]){
	FILE *fp = fopen(filename, "w");
	int i;
	for (i = 1; i <= N; i++){
		
		fprintf(fp, "%Le   ", Pi[i]);
	}
	fclose(fp);
}

//find the average of models
void avgOfAvgModels(int totalIters){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			Aaverage[i][j] /= totalIters;

		}
	}
	for (i = 1; i <= N; i++){
		for (j = 1; j <= M; j++){
			Baverage[i][j] /= totalIters;
		}
	}
	for (i = 1; i <= N; i++){
		Piaverage[i] /= totalIters;
	}
}

//save average model to disk
void saveAvgModel(int digit, char* user){
	char Aavgloc[100], BavgLoc[100], PIavgLoc[100], ind[3];

	sprintf(Aavgloc, "output/savedmodels/%s/digit_%d_A.txt", user, digit);
	FILE *fp = fopen(Aavgloc, "w");
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fprintf(fp, "%Le   ", Aaverage[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	
	sprintf(BavgLoc, "output/savedmodels/%s/digit_%d_B.txt", user, digit);
	fp = fopen(BavgLoc, "w");
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fprintf(fp, "%Le   ", Baverage[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	sprintf(PIavgLoc, "output/savedmodels/%s/digit_%d_PI.txt", user, digit);
	fp = fopen(PIavgLoc, "w");
	for(int i=1; i<=N; i++){
		fprintf(fp, "%Le   ", Piaverage[i]);
	}
	fclose(fp);
}

//preprocess a speech file
void preprocess(char *filename){
    FILE *fp;
    long int totalSample = 0;
    char line[100];

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("Error opening file\n");
    }

	dcshift = 0;
    absmax = 0;
    while(!feof(fp)){
        fgets(line, 100, fp);
        if(!isalpha(line[0])){
            totalSample++;
			dcshift += atoi(line);
            if(absmax < abs(atoi(line)))
                absmax = abs(atoi(line));
        }
    }
    
    multfact = (double)LIMIT/absmax;
    dcshift /= totalSample;

    fclose(fp);
}

//calculating the tokhura distances of a Cis and codebook rows
void calcTokhuraDistance(long double cepstralCoeff[12], int index, FILE *fp){
	int  min_index = 0;
	long double min = DBL_MAX;
	long double sum[CODEBOOKSIZE] = { 0 };
	

	for (int j = 0; j < CODEBOOKSIZE; j++){
		for (int i = 0; i < P; i++){
			sum[j] += tokhuraWeights[i] * (cepstralCoeff[i] - codeBook[j][i])*(cepstralCoeff[i] - codeBook[j][i]);
		}
		if (sum[j] < min){
			min = sum[j];
			min_index = j;
		}
	}

	O[index] = min_index + 1;

	
	fprintf(fp, "%4d ", O[index]);
}

//Calculating Cis
void findCis(){
	
	double sum=0;
	Ci[0]=log(Ri[0]*Ri[0]);

	for(int m=1;m<=P;m++){
		sum=0;
		for(int k=1;k<m;k++){
			sum += (k*Ci[k]*Ai[m-k])/(m*1.0);
		}
		Ci[m]=Ai[m]+sum;
	}
	
}

//Durbin's algorithm to calculate Ais
void durbinAlgo(){
	double alpha[13][13],E[13],K[13];
	double sum = 0;
	E[0] = Ri[0];
	
	for(int i=1;i<=P;i++){
		sum=0;
		for(int j=1;j<=i-1;j++){
			sum += alpha[i-1][j]*Ri[i-j];	
		}
		
		K[i]=(Ri[i]-sum)/E[i-1];
			
		alpha[i][i]=K[i];
	
		for(int j=1;j<=i-1;j++){
			alpha[i][j]=alpha[i-1][j] - K[i]*alpha[i-1][i-j];
		}
	
		E[i]=(1-(K[i]*K[i]))*E[i-1];
	}

	
	for(int i=1;i<=P;i++){
		Ai[i]= alpha[P][i];
	}	
}

//Calculating Ris
void calculate_Ris(double *samp){
	
	for(int m =0; m<=P; m++){
		Ri[m] = 0;
		for(int k=0; k<FRAME_SIZE-m; k++){
			Ri[m] += samp[k]*samp[k+m];
		}
	}
}

//Function to apply raised sin window to Cis
void applyingRaisedSinWindow(){
	long double sum=0;
	for(int m=1;m<=P;m++){
		sum = (P/2)*sin((PI*m)/P);
		Ci[m]*=sum;	
	}	
}

//Calling in order the functions to calculaing the cepstral cofficients
void calcCPrime(double *samp){
	calculate_Ris(samp);
	
	durbinAlgo();
	
	findCis();
	
	applyingRaisedSinWindow();
}

//create the input data
void createInput(){
	
	samplesize = 0;
	
	for(int i=0; i<maxsamplesize; i++){
		inputsamples[samplesize++] = sample[i];
		
	}
	
}

//Generate the observation sequence
void generateObsSeq(char *filename){
	int obs_ind = 1;
	FILE *op = fopen(filename, "w");
	if(op == NULL) {
		printf("Couldn't open file!\n");
		exit(1);
	}
	
	createInput();
	double fsamp[FRAME_SIZE];
	int num_frames = 0;
	for(int i=0; i<samplesize; i++){
		num_frames++;
		for(int j = 0; j<320; j++)
			fsamp[j] = inputsamples[i++];

		calcCPrime(fsamp);
		calcTokhuraDistance(Ci, obs_ind++, op);
	}
	T = num_frames;
	
	fprintf(op, "\n");
	fclose(op);
	
}

//Train our model
void trainModel(char* user){
	char filename[100], line[100], obs_file[100], dump_file[100], com_dump[100];
	deleteModels();
	int total_files_trained = 2;

	for(int d = 0; d<=13; d++){
		printf("Training model for digit %d......................\n",d);
		deleteTheModel();
		for(int u = 1; u <=  total_files_trained; u++){
			sprintf(filename, "dataset/%s/%s_E_%d_%d.txt",user, user, d, u);
			FILE *f = fopen(filename, "r");
			if(f == NULL){
				printf("Issue in opening file %s", filename);
				exit(1);
			}
			preprocess(filename);
			maxsamplesize = 0;
			while(!feof(f)){
				fgets(line, 100, f);
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcshift)*multfact);
					
					sample[maxsamplesize++] = normalizedX;
				}
			}
			fclose(f);
			sprintf(obs_file, "output/obsSeqTrain/HMM_OBS_SEQ_%d_%d.txt", d, u);
			generateObsSeq(obs_file);
			if(train == 0)
				initializeModel(d, 1, "skeletalModel");
			else
				initializeModel(d, 1, "avg");
			int iteration = 1;
			probgivenOandLambda = 0, XprobgivenOandLambda = -1;
			while(probgivenOandLambda > XprobgivenOandLambda && iteration < 1000){
				iteration++;
				XprobgivenOandLambda = probgivenOandLambda; 
				forwardProcedure();
				backwardProcedure();
				viterbiAlgorithm();
				findXi();
				findGamma();
				recalculateModel();
			}
			addToAvgModel();
		}
		avgOfAvgModels(total_files_trained);
		saveAvgModel(d,user);
		deleteAvgModel();
	}
	train++;
}

//Load codebook from disk
void loadCodeGeneral(){
	FILE* fp = fopen("codebook.txt", "r");												
	if(!fp)
	{
		printf("Couldn't read file!\n");
		exit(1);
	}
	while(!feof(fp)){
		for(int i = 0; i < CODEBOOKSIZE; i++){
			for(int j = 0; j < P; j++)
				fscanf(fp, "%Lf", &codeBook[i][j]);
		}
	}

	fclose(fp);
}

//Load codebook from disk
void loadCode(char* user){
	char filename[100];
	sprintf(filename, "codebook/%s/codebook.txt", user);
	FILE* fp = fopen(filename, "r");												
	if(!fp)
	{
		printf("Couldn't read file!\n");
		exit(1);
	}
	while(!feof(fp)){
		for(int i = 0; i < CODEBOOKSIZE; i++){
			for(int j = 0; j < P; j++)
				fscanf(fp, "%Lf", &codeBook[i][j]);
		}
	}

	fclose(fp);
}

//Predict the digit
int getAnswerOperator(char* user){
	answer = 0;
	maximumprob = 0;
	for(int k = 10; k<=13; k++){
		readAvgModel(k,user);
		forwardProcedurePredict(k);
	}

	return answer;
}

//Predict the digit
int getAnswerDigit(char* user){
	answer = 0;
	maximumprob = 0;
	for(int k = 0; k<=9; k++){
		readAvgModel(k,user);
		forwardProcedurePredict(k);
	}

	return answer;
}

//Digit operator from recording
int predictoperator(char* user){
	char obs_file[100], line[100];
	printf("Recording will start now.\n");
	fflush(stdout);
	system("pause");
	system("cls");
	system("Recording_Module.exe 3 input.wav input_file.txt");
	system("cls");
	initializeModel(0, 0, "skeletalModel");
	FILE *f = fopen("input_file.txt", "r");
	if(f == NULL){
		printf("Issue in opening file input_file.txt");
		exit(1);
	}
	preprocess("input_file.txt");
	maxsamplesize = 0;
	while(!feof(f)){
		fgets(line, 100, f);
		if(!isalpha(line[0])){
			int y = atof(line);
			double normalizedX = floor((y-dcshift)*multfact);
			
			sample[maxsamplesize++] = normalizedX;
		}
	}
	fclose(f);
	generateObsSeq("output/obsSeqTestPredict/obs_seq.txt");
	int predicteddigit = getAnswerOperator(user);
	return predicteddigit;
}

//Digit predicting from recording
int predictdigit(char* user){
	char obs_file[100], line[100];
	printf("Recording will start now.\n");
	fflush(stdout);
	system("pause");
	system("cls");
	system("Recording_Module.exe 3 input.wav input_file.txt");
	system("cls");
	initializeModel(0, 0, "skeletalModel");
	FILE *f = fopen("input_file.txt", "r");
	if(f == NULL){
		printf("Issue in opening file input_file.txt");
		exit(1);
	}
	preprocess("input_file.txt");
	maxsamplesize = 0;
	while(!feof(f)){
		fgets(line, 100, f);
		if(!isalpha(line[0])){
			int y = atof(line);
			double normalizedX = floor((y-dcshift)*multfact);
			
			sample[maxsamplesize++] = normalizedX;
		}
	}
	fclose(f);
	generateObsSeq("output/obsSeqTestPredict/obs_seq.txt");
	int predicteddigit = getAnswerDigit(user);
	return predicteddigit;
}

//Model testing
void testModel(char* user){
	char filename[100], line[100], test_file[100];
	int correctAns = 0, totalAns = 0;
	for(int d = 0; d<=13; d++){
		for(int j = 1; j<3; j++){
			sprintf(filename, "dataset/%s/%s_E_%d_%d.txt",user, user, d, j);
			FILE *f = fopen(filename, "r");
			if(f == NULL){
				printf("Issue in opening file input_file.txt");
				exit(1);
			}
			preprocess(filename);
			maxsamplesize = 0;
			while(!feof(f)){
				fgets(line, 100, f);
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcshift)*multfact);
					sample[maxsamplesize++] = normalizedX;
				}
			}
			fclose(f);
			generateObsSeq("output/obsSeqTestPredict/obs_seq.txt");
			answer = 0;
			maximumprob = 0;
			for(int k = 0; k<=13; k++){
				readAvgModel(k, user);
				forwardProcedurePredict(k);
				deleteAvgModel();
			}
			printf("Predicted digit: %d, True Digit: %d\n", answer, d);
			if(answer == d) 
			{
				correctAns++; 
			}
			totalAns++;
		}
	}
	printf("Accuracy: %f %% \n", (correctAns*100.0)/totalAns);
}

int checkuser(char* user)
{
	char filename[100];
	FILE* fp;
	for(int i=0;i<=13;i++)
	{
		sprintf(filename, "output/savedmodels/%s/digit_%d_A.txt", user,i);
		fp=fopen(filename,"r");
		if(fp==NULL)
		{
			return 0;
		}
		fclose(fp);

		sprintf(filename, "output/savedmodels/%s/digit_%d_B.txt", user,i);
		fp=fopen(filename,"r");
		if(fp==NULL)
		{
			return 0;
		}
		fclose(fp);

		sprintf(filename, "output/savedmodels/%s/digit_%d_PI.txt", user,i);
		fp=fopen(filename,"r");
		if(fp==NULL)
		{
			return 0;
		}
		fclose(fp);
	}
	return 1;
}

int checkuserdata(char* user)
{
	char filename[100];
	FILE* fp;
	for(int i=0;i<=13;i++)
	{
		for(int j=1;j<3;j++)
		{
			sprintf(filename, "dataset/%s/%s_E_%d_%d.txt", user, user, i,j);
			fp=fopen(filename,"r");
			if(fp==NULL)
			{
				return 0;
			}
			fclose(fp);
		}
		
	}
	return 1;
}

void recorduserdata(char* user)
{
	
	char command[100];
	sprintf(command, "dataset\\%s",user);
	_mkdir(command);
	sprintf(command, "output\\savedmodels\\%s",user);
	_mkdir(command);
	for(int i=0;i<10;i++)
	{
		for(int j=1;j<3;j++)
		{
			printf("Utter the %dth utterence for the digit %d\n",j,i);
			printf("Press any key to start recoding.(the recording will go for 2 seconds.)\n");
			fflush(stdout);
			system("pause");
			system("cls");
			sprintf(command, "Recording_Module.exe 2 input.wav dataset/%s/%s_E_%d_%d.txt", user,user,i,j);
			system(command);
			system("cls");
		}
	}

	for(int j=1;j<3;j++)
	{
		printf("Utter the %dth utterence for the operation +\n",j);
		printf("Press any key to start recoding.(the recording will go for 2 seconds.)\n");
		fflush(stdout);
		system("pause");
		system("cls");
		sprintf(command, "Recording_Module.exe 2 input.wav dataset/%s/%s_E_10_%d.txt", user,user,j);
		system(command);
		system("cls");
	}

	for(int j=1;j<3;j++)
	{
		printf("Utter the %dth utterence for the operation -\n",j);
		printf("Press any key to start recoding.(the recording will go for 2 seconds.)\n");
		fflush(stdout);
		system("pause");
		system("cls");
		sprintf(command, "Recording_Module.exe 2 input.wav dataset/%s/%s_E_11_%d.txt", user,user,j);
		system(command);
		system("cls");
	}

	for(int j=1;j<3;j++)
	{
		printf("Utter the %dth utterence for the operation *\n",j);
		printf("Press any key to start recoding.(the recording will go for 2 seconds.)\n");
		fflush(stdout);
		system("pause");
		system("cls");
		sprintf(command, "Recording_Module.exe 2 input.wav dataset/%s/%s_E_12_%d.txt", user,user,j);
		system(command);
		system("cls");
	}

	for(int j=1;j<3;j++)
	{
		printf("Utter the %dth utterence for the operation \\\n",j);
		printf("Press any key to start recoding.(the recording will go for 2 seconds.)\n");
		fflush(stdout);
		system("pause");
		system("cls");
		sprintf(command, "Recording_Module.exe 2 input.wav dataset/%s/%s_E_13_%d.txt", user,user,j);
		system(command);
		system("cls");
	}
}

int _tmain(int argc, _TCHAR* argv[]){
	loadCodeGeneral();
	while(true){
		int choice,choice2;
		char username[100];
		int useravail, operand1, operand2, operatorsym; 
		printf("Enter your choice:\n");
		printf("1.\tCalculate\n");
		printf("2.\tTrain registered user's model again\n");
		printf("3.\tExit\n");
		fflush(stdout);
		scanf("%d",&choice);
		system("cls");
		if(choice==1)
		{
			printf("Please Enter you username.\n");
			fflush(stdout);
			scanf("%s",username);
			system("cls");
	 		useravail = checkuser(username);
			if(useravail)
			{
				printf("Hi %s!\n",username);
				printf("Utter the first operand (0-9):\n");
				operand1 = predictdigit(username);
				system("cls");
				printf("Utter the operator (+,-,*,/):\n");
				operatorsym = predictoperator(username);
				system("cls");
				printf("Utter the second operand (0-9):\n");
				operand2 = predictdigit(username);
				system("cls");
				if(operatorsym==10)
				{
					printf("%d + %d = %d\n",operand1,operand2,operand1+operand2);
					FILE *fg=fopen("output.txt","w");
					fprintf(fg,"%d + %d = %d",operand1,operand2,operand1+operand2);
					fclose(fg);
					system("python speak.py");
				}
				else if(operatorsym==11)
				{
					printf("%d - %d = %d\n",operand1,operand2,operand1-operand2);
					FILE *fg=fopen("output.txt","w");
					fprintf(fg,"%d minus %d = %d",operand1,operand2,operand1-operand2);
					fclose(fg);
					system("python speak.py");
				}
				else if(operatorsym==12)
				{
					printf("%d * %d = %d\n",operand1,operand2,operand1*operand2);
					FILE *fg=fopen("output.txt","w");
					fprintf(fg,"%d multiply %d = %d",operand1,operand2,operand1*operand2);
					fclose(fg);
					system("python speak.py");
				}
				else if(operatorsym==13)
				{
					if(operand2==0)
					{
			 			printf("%d / %d = ? Cannot divide by zero\n",operand1,operand2);
						FILE *fg=fopen("output.txt","w");
						fprintf(fg,"Cannot divide by zero");
						fclose(fg);
						system("python speak.py");
					}
			 		else
					{
			 			printf("%d / %d = %f\n",operand1,operand2,float(operand1)/float(operand2));
						FILE *fg=fopen("output.txt","w");
						fprintf(fg,"%d divide by %d = %f\n",operand1,operand2,float(operand1)/float(operand2));
						fclose(fg);
						system("python speak.py");
					}
				}
				system("pause");
				system("cls");
			}
			else{
				int isdataready=checkuserdata(username);
				if(isdataready){
					printf("%s, Your dataset is available but trained model is not available. Should we collect your recording or train your model from the present dataset?\n",username);
					printf("Enter your choice:\n");
					printf("1.\tCollect my recordings.\n");
					printf("2.\tTrain already available data.\n");
					printf("3.\tGo to main menu\n");
					fflush(stdout);
					scanf("%d",&choice2);
					system("cls");
					if(choice2==2)
					{
						trainModel(username);
						testModel(username);
						system("pause");
					}
					else if(choice2==1){
						recorduserdata(username);
						trainModel(username);
						testModel(username);
						system("pause");
					}
					
				}
				else{
					printf("User %s is not registered!\n",username);
					printf("Register User?\n");
					printf("Enter your choice:\n");
					printf("1.\tRegister User\n");
					printf("2.\tGo to main menu\n");
					fflush(stdout);
					scanf("%d",&choice2);
					system("cls");
					if(choice2==1)
					{
						recorduserdata(username);
						trainModel(username);
						testModel(username);
						system("pause");
					}
				}
				
			}
		}
		else if(choice==2){
			printf("Please Enter you username.\n");
			fflush(stdout);
			scanf("%s",username);
			system("cls");
			int isdataready=checkuserdata(username);
				if(isdataready){
					printf("%s, Your dataset is available but trained model is not available. Should we collect your recording or train your model from the present dataset?\n",username);
					printf("Enter your choice:\n");
					printf("1.\tCollect my recordings.\n");
					printf("2.\tTrain already available data.\n");
					printf("3.\tGo to main menu\n");
					fflush(stdout);
					scanf("%d",&choice2);
					system("cls");
					if(choice2==2)
					{
						trainModel(username);
						testModel(username);
						system("pause");
					}
					else if(choice2==1){
						recorduserdata(username);
						trainModel(username);
						testModel(username);
						system("pause");
					}
					
				}
				else{
					
					
						recorduserdata(username);
						trainModel(username);
						testModel(username);
						system("pause");
					
				}
		}
		else{
			printf("Thank you!");
	 		exit(1);
		}
		system("cls");
		
		
	}
	
	return 0;
}