#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <functional>
#include <utility>
#include <math.h>

//This file contains a simulation of an SSN network with 4 cell types. The connection strengths used are found by using least squares approach and the size tuning curves of the inputs from L4 and LM are taken from data. We want to see if the fixed point of the system (response profiles) is consistent with the data and play with counterfactuals.

#define PBC 1//periodic boundary conditions, 0 for open boundaries

using namespace std;

const int N=30; //N*N units of each cell type.
const double dx2deg2=6*6;//(degrees per lattice spacing) squared
//N*sqrt(dx2deg2)=tot visual field considered (in degrees)

struct fr{
	double *Pyr = new double[N*N];//Pyr firing rates
	double *PV = new double[N*N];
	double *SOM = new double[N*N];
	double *VIP = new double[N*N];
} FiringRatesL3;


struct paramt{ //parameters of the inputs, as found by the fits or parametrizations (see file fit_sizetunings.m)
	double bd;
	double kdc;
	double kd;
	double c1d;
	double c2d;
	double s1d;
	double s2d;

	double bi;
	double kic1;
	double kic2;
	double ki;
	double c1i1;
	double c1i2;
	double c2i;
	double s1i;
	double s2i;

} Prmt4,PrmtM;


struct jij{ //connectivity matrix
	double PyrPyr;
	double PyrPV;
	double PyrSOM;
	double L4toPyr;
	double LMtoPyr;

	double PVPyr;
	double PVPV;
	double PVSOM;
	double L4toPV;
	double LMtoPV;

	double SOMPyr;
	double SOMPV;
	double SOMVIP;
	double LMtoSOM;
	double LMextratoSOM;

	double VIPPyr;
	double VIPPV;
	double VIPSOM;
	double LMtoVIP;

	//double SOMSOM;//double PyrVIP;//double VIPVIP;//double PVVIP; // these are deliberately hardcoded to zero.

	double sigma2PyrPyr2=2*64;//theese are 2sigma_{AB}^2
	double sigma2PVPV2=2*25;//
	double sigma2PyrPV2=2*25;//
	double sigma2PyrSOM2=2*25;
	double sigma2PVPyr2=2*64;
	double sigma2PVSOM2=2*25;
	double sigma2SOMPyr2=2*100;
	double sigma2SOMPV2=2*25;
	double sigma2SOMVIP2=2*25;//
	double sigma2VIPPyr2=2*64;
	double sigma2VIPPV2=2*25;
	double sigma2VIPSOM2=2*25;

	double biasPyr;
	double biasPV;
	double biasSOM;
	double biasVIP;
} J;

double *H4 = new double[N*N];
double *HLM = new double[N*N];
double *HLMextra = new double[N*N];

double **WPyrPyr = new double*[N*N];
double **WPVPV = new double*[N*N];	
double **WPyrPV = new double*[N*N];//WPyrPV[i][j] is the strength of the connection from PV neuron j to Pyr neuron i
double **WPyrSOM = new double*[N*N];

double **WPVPyr = new double*[N*N];
double **WPVSOM = new double*[N*N];
double **WSOMPyr = new double*[N*N];
double **WSOMPV = new double*[N*N];
double **WSOMVIP = new double*[N*N];

double **WVIPPyr = new double*[N*N];
double **WVIPPV = new double*[N*N];
double **WVIPSOM = new double*[N*N];

//SOM AND VIP ARE SLOWER AND DO NOT RECEIVE INPUT FROM L4 BUT DO RECEIVE INPUT FROM LM
double tauPyr=0.03;//this is 30ms--> the rates will be in Hertz
double tauPV=0.01;//10ms, Fast Spiking cells
double tauSOM=0.02;//20ms
double tauVIP=0.02;//20ms
double k = 1;//constant for the transfer function (can be reabsorbed in the WAB and biases)
double n = 2;//exponent of the power law nonlinearity in the transfer function

void InitializeConnectionsAndBias(jij &);
void ReadParamL4(paramt &);
void ReadParamLM(paramt &);
void InitializeL3(void);
void DefineL4Inputs(int,int);
void DefineLMInputs(int,int);
void Integrate2DPowerLawIO(double);
double finddistance(int,int);
void CompareWithAnalytics(void);

int center=round(N*(N-1)/2);
int top=round(N/2);

int Htype;
//HTYPE SPECIFIES THE STIMULUS CONDITION. 1 for direct 5 for inverse


int main(int argc, char *argv[])
{

	// ################ PARAMETERS INTRODUCED FROM TERMINAL ################ 
	int T=atoi(argv[1]);//number of simulation steps
	int stimsize=atoi(argv[2]);//stimulus size in degrees (used for the size tuning curves of LM and L4 inputs)
	Htype=atoi(argv[3]);//SPECIFIES THE STIMULUS CONDITION. 1 for direct 5 for inverse
	int thread=atoi(argv[4]);//just to distinguish datasets with the same parameters

	srand (time(NULL));
	int tstep=100;

	// ################ INTEGRATION VARIABLES ################ 
	double dt=0.001;
	int i,j,t,nn;
	double d2,d;

	time_t start,end;
	time (&start);


	// ################ OUTPUT FILES ################ 
	char Name1[100],Name2[100],Name3[100],Name4[100],Name5[100],Name6[100],Name7[100],Name8[100],Name9[100],Name10[100],Param[100],Name11[100],Name14[100];
	char Name12[100],Name16[100];
	//depending on what you want to analyze you should give more detailed names or save parameter value
	if(Htype==1)
		sprintf(Param,"Dirsiz%d_%d.txt",stimsize,thread); 
	else
		sprintf(Param,"Invsiz%d_%d.txt",stimsize,thread); 

	strcpy(Name1,"FR_ts_");
	strcat(Name1,Param);
	strcpy(Name2,"AllFR_Pyr_");
	strcat(Name2,Param);
	strcpy(Name5,"AllFR_PV_");
	strcat(Name5,Param);
	strcpy(Name6,"AllFR_SOM_");
	strcat(Name6,Param);
	strcpy(Name7,"AllFR_VIP_");
	strcat(Name7,Param);
	strcpy(Name3,"FFInputL4_");
	strcat(Name3,Param);
	strcpy(Name4,"FFinputLM_");
	strcat(Name4,Param);

	fstream RTF1(Name1),RTF2(Name2),RTF3(Name3),RTF4(Name4),RTF5(Name5),RTF6(Name6),RTF7(Name7),RTF10(Name10),RTF12(Name12);
	RTF1.open(Name1, ios::out | ios::trunc); 
	RTF2.open(Name2, ios::out | ios::trunc); 
	RTF3.open(Name3, ios::out | ios::trunc); 
	RTF4.open(Name4, ios::out | ios::trunc);
	RTF5.open(Name5, ios::out | ios::trunc);
	RTF6.open(Name6, ios::out | ios::trunc); 
	RTF7.open(Name7, ios::out | ios::trunc); 

	// ################ ALLOCATE MEMORY FOR THE CONNECTIVITY MATRIX ################ 
	for(int i=0;i<N*N;i++){
		WPyrPyr[i] = new double[N*N];
		WPVPV[i] = new double[N*N];	
		WPyrPV[i] = new double[N*N];//WPyrPV[i][j] is the strength of the connection from PV neuron j to Pyr neuron i
		WPyrSOM[i] = new double[N*N];

		WPVPyr[i] = new double[N*N];
		WPVSOM[i] = new double[N*N];	
		WSOMPyr[i] = new double[N*N];
		WSOMPV[i] = new double[N*N];
		WSOMVIP[i] = new double[N*N];

		WVIPPyr[i] = new double[N*N];
		WVIPPV[i] = new double[N*N];
		WVIPSOM[i] = new double[N*N];
	}

	double averageFRPyr=0,RecInputE=0, RecInputI=0;

	// ################ INITIALIZATIONS ################ 
	InitializeL3();
	cout << "\tNetwork Initialized" << endl;
	InitializeConnectionsAndBias(J);
	cout << "\tInteractions and Biases Initialized" << endl;

	DefineL4Inputs(Htype,stimsize);
	for(i=0;i<N*N;i++) {
		RTF3 << H4[i] << " ";
		if((i+1)%N==0)
			RTF3 << "\n";
	}
	RTF3.close();
	cout << "\tFF input Initialized" << endl;

	DefineLMInputs(Htype,stimsize);
	for(i=0;i<N*N;i++) {
		RTF4 << HLM[i] << " ";
		if((i+1)%N==0)
			RTF4 << "\n";
	}
	RTF4.close();
	cout << "\tFB input Initialized" << endl;

	CompareWithAnalytics();//check these values with the file compare_simulationsanalytics.nb, when the ConM imported is the same!

	// ################ INTEGRATION ################ 
	cout << "\tIntegration Started" << endl;

	for(t=0;t<T;t++){	

		Integrate2DPowerLawIO (dt); //<<integration happens here! LM ACTIVE

		//output firing rate in the center and on the top edge at time t -- TO CHECK CONVERGENCE
		RTF1 << FiringRatesL3.Pyr[center] << " " << FiringRatesL3.PV[center] << " " << FiringRatesL3.SOM[center] << " " << FiringRatesL3.VIP[center] << endl;
		nn=round(100*t/(T-1));
		cout << "\r\t" << nn << "\% complete";
		fflush(stdout);

	}
	// ################ END OF INTEGRATION ################ 
	cout << "\tThe End" << endl;

	cout << "\tWriting spatial profiles" << endl;
	//output firing rates of all L3 units at final time T
	for(i=0;i<N*N;i++){ 
		RTF2 << FiringRatesL3.Pyr[i] << " ";
		RTF5 << FiringRatesL3.PV[i] << " ";
		RTF6 << FiringRatesL3.SOM[i] << " ";
		RTF7 << FiringRatesL3.VIP[i] << " ";
		if((i+1)%N==0) {
			RTF2 << "\n";
			RTF5 << "\n";
			RTF6 << "\n";
			RTF7 << "\n";
		}
	}
	RTF2.close();
	RTF5.close();
	RTF6.close();
	RTF7.close();

	time(&end);

	CompareWithAnalytics();
	cout<<"Elapsed time is " << end-start << " s\n";
	//system("say -v Vicki i am done");
	return EXIT_SUCCESS;

}









// ####################################################################
// ################ DEFINITION OF THE USED FUNCTIONS ################ 
// ####################################################################



//set initial values of firing rates at baseline level
void InitializeL3(void){
	//uniform_real_distribution<double> ran_u(0.0,1.0); //Uniform number distribution between 0 and 1 Hz
	//mt19937 gen(rand());; 
	double r;

	for(int i=0;i<N*N;i++){
		r = 0.5996;//ran_u(gen);
		FiringRatesL3.Pyr[i]=r;
		r = 1.35;//ran_u(gen);
		FiringRatesL3.PV[i]=r;
		r = 0.6818;//ran_u(gen);
		FiringRatesL3.SOM[i]=r;
		r = 0.8437;//ran_u(gen);
		FiringRatesL3.VIP[i]=r;
	}
}



//define the INPUT CURRENT FROM L4 (i.e. a convolution between activity and projections)
void DefineL4Inputs(int Htype, int stimsize){

	//input parameters as defined by the experimental measures (fits/parametrization)
	ReadParamL4(Prmt4);
	double dd;
	double sigmaL4_2=2.*(Prmt4.kdc+Prmt4.kd*stimsize)*(Prmt4.kdc+Prmt4.kd*stimsize);
	double sigmaL4i_2=2.*(Prmt4.kic1+Prmt4.ki*stimsize)*(Prmt4.kic1+Prmt4.ki*stimsize);
	double sigmaL4bari_2=2.*(Prmt4.kic2+Prmt4.ki*stimsize)*(Prmt4.kic2+Prmt4.ki*stimsize);
	double sigmaL4toL23_2=2.*8.*8.;
	J.L4toPyr=J.L4toPyr/(3.141592*sigmaL4toL23_2);//Gaussian connectivity is normalized
	J.L4toPV=J.L4toPV/(3.141592*sigmaL4toL23_2);

	//define spatial dependence of the FeedForward current from layer 4
	double *L4 = new double[N*N];

	if(Htype==1){
		//CLASSICAL INPUT is a 2D gaussian field in L4 (see Supplement)
		double A4=Prmt4.c1d*erf(stimsize/Prmt4.s1d)-Prmt4.c2d*erf(stimsize/Prmt4.s2d);
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			L4[i]=1*A4*exp(-(dd*dd*dx2deg2/sigmaL4_2))+Prmt4.bd;//rate field of L4
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}

		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=L4[j]*exp(-(dd*dd*dx2deg2/sigmaL4toL23_2));
			}
			H4[i]=dx2deg2*sum;//input current from L4
		}	
	}

	if(Htype==2){
		//INVERSE INPUT is a difference of 2 gaussians in L4 (see Supplement)
		double A4=Prmt4.c1i1*(erf(stimsize/Prmt4.s1i)-erf(stimsize/Prmt4.s2i));
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			L4[i]=1*A4*(exp(-(dd*dd*dx2deg2/(sigmaL4i_2)))-exp(-(dd*dd*dx2deg2/(sigmaL4bari_2))))+Prmt4.bi;//rate field of L4
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=L4[j]*exp(-(dd*dd*dx2deg2/(sigmaL4toL23_2)));
			}
			H4[i]=dx2deg2*sum;//input current from L4
		}	
	}

	delete[] L4;
}



//define the INPUT CURRENT FROM LM (i.e. a convolution between activity and projections)
void DefineLMInputs(int Htype, int stimsize){

	ReadParamLM(PrmtM);
	//define spatial dependence of the FeedForward current from LM
	double dd;
	double sigmaLM_2=2.*(PrmtM.kdc+PrmtM.kd*stimsize)*(PrmtM.kdc+PrmtM.kd*stimsize);

	double sigmaLMi_2=2.*(PrmtM.kic1+PrmtM.ki*stimsize)*(PrmtM.kic1+PrmtM.ki*stimsize);
	double sigmaLMbari_2=2.*(PrmtM.kic2+PrmtM.ki*stimsize)*(PrmtM.kic2+PrmtM.ki*stimsize);

	double sigmaLMtoL23_2=2.*20.*20.;

	if(stimsize==0){//this is just to avoid patologies in the Gaussian function
		sigmaLM_2=20;
		sigmaLMi_2=20;
		sigmaLMbari_2=20;
	}

	double *LM = new double[N*N];
	double *LMextra = new double[N*N]; //small excitatory extra input rate field do correct for excitatory input to SOM
	J.LMtoPyr=J.LMtoPyr/(3.141592*sigmaLMtoL23_2);
	J.LMtoPV=J.LMtoPV/(3.141592*sigmaLMtoL23_2);
	J.LMtoSOM=J.LMtoSOM/(3.141592*sigmaLMtoL23_2);
	J.LMextratoSOM=J.LMextratoSOM/(3.141592*2*8*8);
	J.LMtoVIP=J.LMtoVIP/(3.141592*sigmaLMtoL23_2);


	//double extraALM=0;
	//if(stimsize<30)
   	//	extraALM=-1.;//.*min(stimsize,55);
	//else if(stimsize<30)
	//	extraALM=0.02*stimsize;
	//else
	//	extraALM=(stimsize-25)*0.015; 


	if(Htype==1){
		//CLASSICAL INPUT is a 2D gaussian field in L4 (see Supplement)
		double ALM=PrmtM.c1d*erf(stimsize/PrmtM.s1d)-PrmtM.c2d*erf(stimsize/PrmtM.s2d);
		double extraALM=sqrt(stimsize);//4*tanh(.3*stimsize);
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			LM[i]=1*ALM*exp(-(dd*dd*dx2deg2/(sigmaLM_2)))+PrmtM.bd;//rate field of LM
			//CHANGE CONSTANT HERE FOR CONTRAST OR OPTOGENETIC MODULATIONS
			LMextra[i]=extraALM*exp(-(dd*dd*dx2deg2/(2*20*20)));
		}

		for(int i=0;i<N*N;i++){
			double sum=0;
			double sum1=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=LM[j]*exp(-(dd*dd*dx2deg2/(sigmaLMtoL23_2)));
				sum1+=LMextra[j]*exp(-(dd*dd*dx2deg2/(2*8*8)));
			}
			HLM[i]=dx2deg2*sum;//input current from LM
			HLMextra[i]=dx2deg2*sum1;
		}	
	}
	if(Htype==2){
		//INVERSE INPUT is a difference of 2 gaussians in L4 (see Supplement)
		double ALM1=PrmtM.c1i1+PrmtM.c1i2*stimsize;
		double ALM2=ALM1-PrmtM.c2i*(erf(stimsize/PrmtM.s1i)-erf(stimsize/PrmtM.s2i));
		cout<<"ALM1="<<ALM1<<endl;
		cout<<"ALM2="<<ALM2<<endl;
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			LM[i]=1*(ALM1*exp(-(dd*dd*dx2deg2/(sigmaLMi_2)))-ALM2*exp(-(dd*dd*dx2deg2/(sigmaLMbari_2))))+PrmtM.bi;//make sure that sigmabari is smaller
			//CHANGE CONSTANT HERE FOR CONTRAST OR OPTOGENETIC MODULATIONS
		}
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=LM[j]*exp(-(dd*dd*dx2deg2/(sigmaLMtoL23_2)));
			}
			HLM[i]=dx2deg2*sum;//input current from LM
		}	
	}
	cout << "totincurrentLM=" << J.LMtoPyr*HLM[1] << endl;

	delete[] LM;
	delete[] LMextra;
}



//simple Euler algorithm for integration of differential equations
void Integrate2DPowerLawIO(double dt){

	int i,j;
	double TotInputPyr,TotInputPV,TotInputSOM,TotInputVIP,RecInputPyr,RecInputPV,RecInputSOM,RecInputVIP,RPyrss,RPVss,RSOMss,RVIPss;

	//Euler rule//
	//update L23
	for(i=0;i<N*N;i++){
		RecInputPyr=0;
		RecInputPV=0;
		RecInputSOM=0;
		RecInputVIP=0;

		for(j=0;j<N*N;j++){
			RecInputPyr+=WPyrPyr[i][j]*FiringRatesL3.Pyr[j]+WPyrPV[i][j]*FiringRatesL3.PV[j]+WPyrSOM[i][j]*FiringRatesL3.SOM[j];		
			RecInputPV+=WPVPyr[i][j]*FiringRatesL3.Pyr[j]+WPVPV[i][j]*FiringRatesL3.PV[j]+WPVSOM[i][j]*FiringRatesL3.SOM[j];		
			RecInputSOM+=WSOMPyr[i][j]*FiringRatesL3.Pyr[j]+WSOMPV[i][j]*FiringRatesL3.PV[j]+WSOMVIP[i][j]*FiringRatesL3.VIP[j];		
			RecInputVIP+=WVIPPyr[i][j]*FiringRatesL3.Pyr[j]+WVIPPV[i][j]*FiringRatesL3.PV[j]+WVIPSOM[i][j]*FiringRatesL3.SOM[j];		
		}

		TotInputPyr=J.L4toPyr*H4[i]+J.LMtoPyr*HLM[i]+dx2deg2*RecInputPyr+J.biasPyr;
		TotInputPV=J.L4toPV*H4[i]+J.LMtoPV*HLM[i]+dx2deg2*RecInputPV+J.biasPV;
		TotInputSOM=J.LMtoSOM*HLM[i]+J.LMextratoSOM*HLMextra[i]+dx2deg2*RecInputSOM+J.biasSOM;
		//TotInputSOM=J.LMtoSOM*HLM[i]+dx2deg2*RecInputSOM+J.biasSOM;
		TotInputVIP=J.LMtoVIP*HLM[i]+dx2deg2*RecInputVIP+J.biasVIP;

		//rectification
		RPyrss=max(0.,TotInputPyr);
		RPVss=max(0.,TotInputPV);
		RSOMss=max(0.,TotInputSOM);
		RVIPss=max(0.,TotInputVIP);

		//power law nonlinearity
		RPyrss=k*pow(RPyrss,n);
		RPVss=k*pow(RPVss,n);
		RSOMss=k*pow(RSOMss,n);
		RVIPss=k*pow(RVIPss,n);

		FiringRatesL3.Pyr[i]+=dt*(1/tauPyr)*(-FiringRatesL3.Pyr[i]+RPyrss);
		//FiringRatesL3.PV[i]+=dt*(1/tauPV)*(-FiringRatesL3.PV[i]+RPVss);
		FiringRatesL3.PV[i]+=0;
		FiringRatesL3.SOM[i]+=dt*(1/tauSOM)*(-FiringRatesL3.SOM[i]+RSOMss);
		FiringRatesL3.VIP[i]+=dt*(1/tauVIP)*(-FiringRatesL3.VIP[i]+RVIPss);
	}

}



//function to find distance of two points in a 2D lattice of linear size N
double finddistance(int p1,int p2){
	
	int yd, xd, xp1, xp2, yp1, yp2;
	double d;
	xp1=p1%N;
	xp2=p2%N;
	yp1=p1/N;
	yp2=p2/N;

	if(PBC==0){
		xd=abs(xp1-xp2);
		yd=abs(yp1-yp2);
	}
	else{
		xd=min(abs(xp1-xp2),abs(xp1-xp2+N));
		xd=min(xd,abs(xp2-xp1+N));
		yd=min(abs(yp1-yp2),abs(yp1-yp2+N));
		yd=min(yd,abs(yp2-yp1+N));
	}
	d=sqrt(xd*xd+yd*yd);
	return d;
}



void InitializeConnectionsAndBias(jij &J){

	string wrd = "";
	string str, line;
	vector<double> parsedM,parsedT;
	
	//load interaction matrix from result of semianalytical NNLS
	ifstream ifs("inferredWAB14.csv");
	
	while(getline(ifs, line)){
		stringstream ss(line);
		while(getline(ss,str,',')){
			for (auto x : str) 
				wrd=wrd+x;//read next letter

			parsedM.push_back(stod(wrd));//word is complete//convert into number
			wrd="";
		}
	}
	ifs.close();

	/*CHECK THE PARSING
	   k=0;
	   for (int c=0;c<4;c++){
	   for (int i=0;i<7;i++) {
	   cout<<parsedM[k]<<"\t";
	   k++;
	   }
	   cout<<endl;
	   }
	   */
	   

	J.PyrPyr=parsedM[0];
	J.PyrPV=parsedM[1];
	J.PyrSOM=parsedM[2];
	J.L4toPyr=parsedM[4];
	J.LMtoPyr=parsedM[5];
	//6 is LMextra to Pyr=0

	J.PVPyr=parsedM[7];
	J.PVPV=parsedM[8];
	J.PVSOM=parsedM[9];
	J.L4toPV=parsedM[11];
	J.LMtoPV=parsedM[12];
	//13 is LMextra to PV=0


	J.SOMPyr=parsedM[14];
	J.SOMPV=parsedM[15];
	J.SOMVIP=parsedM[17];
	J.LMtoSOM=parsedM[19];
	J.LMextratoSOM=parsedM[20];

	J.VIPPyr=parsedM[21];
	J.VIPPV=parsedM[22];
	J.VIPSOM=parsedM[23];
	J.LMtoVIP=parsedM[26];


	//initializes matrices WXY[i][j] i,j=1..N such that interactions depend on the distance(i,j).
	//a unit interacts with its neighbors (similar RFs).
	int i,j;
	double d, d2;

	//normalizing Jij
	int center=round(N*(N-1)/2);

	for(i=0;i<N*N;i++){
		for(j=0;j<N*N;j++){
			d=finddistance(i,j);
			d2=(double)d*d*dx2deg2;
			WPyrPyr[i][j]=(J.PyrPyr*exp(-d2/J.sigma2PyrPyr2))/(3.141592*J.sigma2PyrPyr2);
			WPyrPV[i][j]=-(J.PyrPV*exp(-d2/J.sigma2PyrPV2))/(3.141592*J.sigma2PyrPV2);//SET TO 0 TO CHECK PV PARADOXICAL
			WPyrSOM[i][j]=-(J.PyrSOM*exp(-d2/J.sigma2PyrSOM2))/(3.141592*J.sigma2PyrSOM2);//MODIFY TO CHECK HYPERPOLARIZATION SOM
			WPVPyr[i][j]=(J.PVPyr*exp(-d2/J.sigma2PVPyr2))/(3.141592*J.sigma2PVPyr2);
			WPVPV[i][j]=-(J.PVPV*exp(-d2/J.sigma2PVPV2))/(3.141592*J.sigma2PVPV2);
			WPVSOM[i][j]=-(J.PVSOM*exp(-d2/J.sigma2PVSOM2))/(3.141592*J.sigma2PVSOM2);
			WSOMPyr[i][j]=(J.SOMPyr*exp(-d2/J.sigma2SOMPyr2))/(3.141592*J.sigma2SOMPyr2);
			WSOMPV[i][j]=-(J.SOMPV*exp(-d2/J.sigma2SOMPV2))/(3.141592*J.sigma2SOMPV2);
			WSOMVIP[i][j]=-(J.SOMVIP*exp(-d2/J.sigma2SOMVIP2))/(3.141592*J.sigma2SOMVIP2);
			WVIPPyr[i][j]=(J.VIPPyr*exp(-d2/J.sigma2VIPPyr2))/(3.141592*J.sigma2VIPPyr2);
			WVIPPV[i][j]=-(J.VIPPV*exp(-d2/J.sigma2VIPPV2))/(3.141592*J.sigma2VIPPV2);
			WVIPSOM[i][j]=-(J.VIPSOM*exp(-d2/J.sigma2VIPSOM2))/(3.141592*J.sigma2VIPSOM2);
		}
	}

	//load biases values from result of semianalytical NNLS
	ifstream ifs2("inferredTA14.csv");
	while(getline(ifs2, line))
		parsedT.push_back(stod(line));
	ifs2.close();

	J.biasPyr=parsedT[0];
	J.biasPV=parsedT[1];
	J.biasSOM=parsedT[2];
	J.biasVIP=parsedT[3];

}



//debugging appendix
void CompareWithAnalytics(void){

	double RecInputPyrPyr,RecInputPyrPV, RecInputPyrSOM, HLMavg, HLMextraavg;
	double RecInputSOMPyr,RecInputSOMPV, RecInputSOMVIP, RecInputSOM,  RecInputSOMPyroff,RecInputSOMPVoff, RecInputSOMVIPoff, RecInputSOMoff;

	RecInputPyrPyr=0.;
	RecInputPyrPV=0.;
	RecInputPyrSOM=0.;
	RecInputSOMPyr=0.;
	RecInputSOMPV=0.;
	RecInputSOMVIP=0.;
	RecInputSOM=0.;
	RecInputSOMPyroff=0.;
	RecInputSOMPVoff=0.;
	RecInputSOMVIPoff=0.;
	RecInputSOMoff=0.;

	for(int j=0;j<N*N;j++){
		HLMextraavg+=HLMextra[j];
		HLMavg+=HLM[j];
		RecInputPyrPyr+=WPyrPyr[center][j]*FiringRatesL3.Pyr[j];
		RecInputPyrPV+=WPyrPV[center][j]*FiringRatesL3.PV[j];
		RecInputPyrSOM+=WPyrSOM[center][j]*FiringRatesL3.SOM[j];
		RecInputSOM+=WSOMPyr[center][j]*FiringRatesL3.Pyr[j]+WSOMPV[center][j]*FiringRatesL3.PV[j]+WSOMVIP[center][j]*FiringRatesL3.VIP[j];
		RecInputSOMoff+=WSOMPyr[1][j]*FiringRatesL3.Pyr[j]+WSOMPV[1][j]*FiringRatesL3.PV[j]+WSOMVIP[1][j]*FiringRatesL3.VIP[j];
		RecInputSOMPyr+=WSOMPyr[center][j]*FiringRatesL3.Pyr[j];
		RecInputSOMPV+=WSOMPV[center][j]*FiringRatesL3.PV[j];
		RecInputSOMVIP+=WSOMVIP[center][j]*FiringRatesL3.VIP[j];
		RecInputSOMPyroff+=WSOMPyr[1][j]*FiringRatesL3.Pyr[j];
		RecInputSOMPVoff+=WSOMPV[1][j]*FiringRatesL3.PV[j];
		RecInputSOMVIPoff+=WSOMVIP[1][j]*FiringRatesL3.VIP[j];

	}
	HLMextraavg=HLMextraavg/(N*N);
	HLMavg=HLMavg/(N*N);


	cout << "PYR: Pyrinput=" << dx2deg2*RecInputPyrPyr << "\tPVinput=" << dx2deg2*RecInputPyrPV << "\tSOMinput=" << dx2deg2*RecInputPyrSOM << "\tL4input=" << J.L4toPyr*H4[center] << "\tLMinput=" <<  J.LMtoPyr*HLM[center] <<  "\tbias="<< J.biasPyr << endl;
	cout << "SOMcen: Pyrinput=" << dx2deg2*RecInputSOMPyr << "\tPVinput=" << dx2deg2*RecInputSOMPV << "\tVIPinput=" << dx2deg2*RecInputSOMVIP << "\taverageLMinput=" <<  J.LMtoSOM*HLMavg << "\tbias="<< J.biasSOM << "\taverageextrainput=" << J.LMtoSOM*HLMextraavg << endl;
	//cout << "SOMoff: Pyrinput=" << dx2deg2*RecInputSOMPyroff << "\tPVinput=" << dx2deg2*RecInputSOMPVoff << "\tVIPinput=" << dx2deg2*RecInputSOMVIPoff << "\tLMinput=" <<  J.LMtoSOM*HLM[1] << "\tbias="<< J.biasSOM<< "\tmax extra bias" << 0.0006*85 << "\tTOTinput=" << dx2deg2*RecInputSOMoff+J.LMtoSOM*HLM[1]+J.biasSOM << endl;

}



//read parameters for L4 rate field from files (see fit_sizetunings.m)
void ReadParamL4(paramt &Prmt4){
	string line;

	if(Htype==1){
		ifstream ifs("parametersL4_direct+.csv");

		ifs>>Prmt4.bd;
		ifs>>Prmt4.kdc;
		ifs>>Prmt4.kd;
		ifs>>Prmt4.c1d;
		ifs>>Prmt4.c2d;
		ifs>>Prmt4.s1d;
		ifs>>Prmt4.s2d;

		ifs.close();
	}else{
		ifstream ifs("parametersL4_inverse+.csv");

		ifs>>Prmt4.bi;
		ifs>>Prmt4.kic1;
		ifs>>Prmt4.kic2;
		ifs>>Prmt4.ki;
		ifs>>Prmt4.c1i1;
		ifs>>Prmt4.c2i;
		ifs>>Prmt4.s1i;
		ifs>>Prmt4.s2i;

		ifs.close();
		cout<<Prmt4.bi<<" "<<Prmt4.kic1<<" "<<Prmt4.kic2<<" "<<Prmt4.ki<<endl<<Prmt4.c1i1<<" "<<Prmt4.c2i<<" "<<Prmt4.s1i<<" "<<Prmt4.s2i<<endl;
	}



}



//read parameters for LM rate field from files (see fit_sizetunings.m)
void ReadParamLM(paramt &PrmtM){
	string line;

	if(Htype==1){
		ifstream ifs("parametersLM_direct+.csv");

		ifs>>PrmtM.bd;
		ifs>>PrmtM.kdc;
		ifs>>PrmtM.kd;
		ifs>>PrmtM.c1d;
		ifs>>PrmtM.c2d;
		ifs>>PrmtM.s1d;
		ifs>>PrmtM.s2d;

		ifs.close();

	}else{
		ifstream ifs("parametersLM_inverse+.csv");

		ifs>>PrmtM.bi;
		ifs>>PrmtM.kic1;
		ifs>>PrmtM.kic2;
		ifs>>PrmtM.ki;
		ifs>>PrmtM.c1i1;
		ifs>>PrmtM.c1i2;
		ifs>>PrmtM.c2i;
		ifs>>PrmtM.s1i;
		ifs>>PrmtM.s2i;

		//PrmtM.ki=0.;//FOR COUNTERFACTUAL MODIFICATION inverse LM WIDTH
		//PrmtM.kd=0.;//FOR COUNTERFACTUAL MODIFICATION direct LM WIDTH

		ifs.close();
		cout<<PrmtM.bi<<" "<<PrmtM.kic1<<" "<<PrmtM.kic2<<" "<<PrmtM.ki<<endl<<PrmtM.c1i1<<" "<<PrmtM.c1i2<<" "<<PrmtM.c2i<<" "<<PrmtM.s1i<<" "<<PrmtM.s2i<<endl;
	}

}

