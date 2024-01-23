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

#define PBC 1//periodic boundary conditions

using namespace std;

const int N=60; //N*N units of each cell type.
const double dx2deg2=6*6;
//N*sqrt(dx2deg2)=tot visual field considered (in degrees)

struct fr{
	double *Pyr = new double[N*N];
} FiringRatesL3;


struct jij{
	double PyrPyr=-.4;
	double L4toPyr=1.;
	double LMtoPyr=.75;
	double SOMtoPyr=-.15;
	double sigma2PyrPyr2=2.*49;//theese are 2sigma^2
	double sigmaL4toL23_2=2.*49;
	double sigmaLMtoL23_2=2.*225;
	double sigmaSOMtoPyr_2=2.*16;

} J;

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
	double A1;
	double A2;
	double sa1;
	double sa2;
	double rho1;
	double rho2;
	double srho1;
	double srho2;

} Prmt4,PrmtM,PrmtS;


double *H4 = new double[N*N];
double *HLM = new double[N*N];
double *HS = new double[N*N];

double **WPyrPyr = new double*[N*N];

//SOM AND VIP ARE SLOWER AND DO NOT RECEIVE INPUT FROM L4 BUT DO RECEIVE INPUT FROM LM
double tauPyr=0.01;//this is 10ms--> the rates will be in Hertz
double k = 1;
double n = 2;//nonlinearity


void InitializeConnections(jij &);
void InitializeL3(void);
void DefineL4Inputs(int,int);
void DefineLMInputs(int,int);
void DefineSOMInputs(int,int);
void ReadParamL4(paramt &);
void ReadParamLM(paramt &);
void ReadParamSOM(paramt &);
void Integrate2DPowerLawIO(double);
double finddistance(int,int);
void CompareWithAnalytics(void);

int center=round(N*(N-1)/2);
int top=round(N/2);
int Htype;

int main(int argc, char *argv[])
{
	
	// ################ PARAMETERS INTRODUCED FROM TERMINAL ################ 
	int T=atoi(argv[1]);//number of simulation steps
	int stimsize=atoi(argv[2]);//stimulus size in degrees (used for the size tuning curves of LM and L4 inputs)
	Htype=atoi(argv[3]);//SPECIFIES THE STIMULUS CONDITION. 1 for direct 5 for inverse
	int thread=atoi(argv[4]);//just to distinguish datasets with the same parameters

	srand (time(NULL));
	int tstep=100;

	//INTEGRATION VARIABLES
	double dt=0.001;
	int i,j,t,nn;
	double d2,d;

	time_t start,end;
	time (&start);


	char Name1[100],Name2[100],Name3[100],Name4[100],Name5[100],Name6[100],Name7[100],Name8[100],Name9[100],Name10[100],Param[100],Name11[100],Name14[100];
	char Name12[100],Name16[100];
	//depending on what you want to analyze you should give more detailed names or save parameter value
	
	if(Htype==1)
		sprintf(Param,"Dirsiz%d_%d.txt",stimsize,thread);  
	else if(Htype==2)
		sprintf(Param,"Invsiz%d_%d.txt",stimsize,thread); 


	strcpy(Name1,"./1ct/FR_ts_");
	strcat(Name1,Param);
	strcpy(Name2,"./1ct/AllFR_Pyr_");
	strcat(Name2,Param);
	strcpy(Name3,"./1ct/FFInputL4_");
	strcat(Name3,Param);
	strcpy(Name4,"./1ct/FFInputLM_");
	strcat(Name4,Param);
	strcpy(Name5,"./1ct/FFInputSOM_");
	strcat(Name5,Param);



	fstream RTF1(Name1),RTF2(Name2),RTF3(Name3),RTF4(Name4),RTF5(Name5),RTF6(Name6),RTF7(Name7),RTF10(Name10),RTF12(Name12);
	RTF1.open(Name1, ios::out | ios::trunc); 
	RTF2.open(Name2, ios::out | ios::trunc); 
	RTF3.open(Name3, ios::out | ios::trunc); 
	RTF4.open(Name4, ios::out | ios::trunc); 
	RTF5.open(Name5, ios::out | ios::trunc); 

	for(int i=0;i<N*N;i++)
		WPyrPyr[i] = new double[N*N];

	double averageFRPyr=0,RecInputE=0, RecInputI=0;

	//>>>>>>>>>>>>>>>>>>>>>>>>INITIALIZATIONS
	InitializeL3();
	cout << "\tNetwork Initialized" << endl;
	InitializeConnections(J);
	cout << "\tInteractions Initialized" << endl;

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

	DefineSOMInputs(Htype,stimsize);
	for(i=0;i<N*N;i++) {
		RTF5 << HS[i] << " ";
		if((i+1)%N==0)
			RTF5 << "\n";
	}
	RTF5.close();
	cout << "\tSOM input Initialized" << endl;


	//if(stimsize==0) 
	CompareWithAnalytics();//check these values with the section debugging in the file AllFiguresCode.nb
	
	//>>>>>>>>>>>>>>>>>>>>>>>>INTEGRATION
	cout << "\tIntegrating" << endl;

	for(t=0;t<T;t++){	
		Integrate2DPowerLawIO (dt); //<<integration happens here!
			
		//output firing rate in the center and on the top edge at time t
		RTF1 << FiringRatesL3.Pyr[center] << endl;
		nn=round(100*t/(T-1));
		cout << "\r\t" << nn << "\% complete";
		fflush(stdout);

	}
	//>>>>>>>>>>>>>>>>>>>>>>>>END OF INTEGRATION
	cout << "\tThe End" << endl;

	cout << "\tWriting spatial profiles" << endl;

	//output firing rates of all L3 units at final time T
	for(i=0;i<N*N;i++){ 
		RTF2 << FiringRatesL3.Pyr[i] << " ";
		if((i+1)%N==0) 
			RTF2 << "\n";
	}
	RTF2.close();

	time(&end);
	CompareWithAnalytics();//check these values with the section debugging in the file AllFiguresCode.nb

	cout<<"Elapsed time is " << end-start << " s\n";
	//system("say -v Vicki i am done");
	return EXIT_SUCCESS;

}



void InitializeL3(void){

	//uniform_real_distribution<double> ran_u(0.0,0.1); //Uniform number distribution
	//mt19937 gen(rand());; 
	double r;

	r = 0.5996;//ran_u(gen);
	for(int i=0;i<N*N;i++)
		FiringRatesL3.Pyr[i]=r;
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
		ifs>>Prmt4.A1;
		ifs>>Prmt4.rho1;
		ifs>>Prmt4.srho1;
		ifs>>Prmt4.srho2;

		ifs.close();
		//cout<<Prmt4.kic1<<" "<<Prmt4.kic2<<" "<<Prmt4.ki<<endl<<Prmt4.A1<<" "<<Prmt4.rho1<<" "<<Prmt4.srho1<<" "<<Prmt4.srho2<<endl;
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
		ifs>>PrmtM.A1;
		ifs>>PrmtM.A2;
		ifs>>PrmtM.rho1;
		ifs>>PrmtM.srho1;
		ifs>>PrmtM.srho2;

		//PrmtM.ki=0.;//FOR COUNTERFACTUAL MODIFICATION inverse LM WIDTH
		//PrmtM.kd=0.;//FOR COUNTERFACTUAL MODIFICATION direct LM WIDTH

		ifs.close();
		//cout<<" "<<PrmtM.kic1<<" "<<PrmtM.kic2<<" "<<PrmtM.ki<<endl<<PrmtM.A1<<" "<<PrmtM.A2<<" "<<PrmtM.rho1<<" "<<PrmtM.srho1<<" "<<PrmtM.srho2<<endl;
	}

}



//read parameters for L4 rate field from files (see fit_sizetunings.m)
void ReadParamSOM(paramt &PrmtS){
	string line;

	if(Htype==1){
		ifstream ifs("parametersSOM_direct+.csv");
		ifs>>PrmtS.bd;
		ifs>>PrmtS.kdc;
		ifs>>PrmtS.kd;
		ifs>>PrmtS.c1d;
		ifs>>PrmtS.c2d;
		ifs>>PrmtS.s1d;
		ifs>>PrmtS.s2d;

		ifs.close();
	}else{

		ifstream ifs("parametersSOM_inverse+.csv");
		ifs>>PrmtS.bi;
		ifs>>PrmtS.kic1;
		ifs>>PrmtS.kic2;
		ifs>>PrmtS.ki;
		ifs>>PrmtS.A1;
		ifs>>PrmtS.A2;
		ifs>>PrmtS.sa1;
		ifs>>PrmtS.sa2;
		ifs>>PrmtS.rho1;
		ifs>>PrmtS.rho2;
		ifs>>PrmtS.srho1;
		ifs>>PrmtS.srho2;

		ifs.close();
		//cout<<PrmtS.kic1<<" "<<PrmtS.kic2<<" "<<PrmtS.ki<<endl<<PrmtS.A1S<<" "<<PrmtS.A2S<<" "<<PrmtS.sa1S<<" "<<PrmtS.sa2S<<" "<<PrmtS.rho1S<<" "<<PrmtS.rho2S<<" "<<PrmtS.srho1S<<" "<<PrmtS.srho2S<<endl;
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
	J.L4toPyr=J.L4toPyr/(3.141592*J.sigmaL4toL23_2);//Gaussian connectivity is normalized

	//define spatial dependence of the FeedForward current from layer 4
	double *L4 = new double[N*N];

	if(Htype==1){
		//CLASSICAL INPUT is a 2D gaussian field in L4 (see Supplement)
		double A4=Prmt4.c1d*erf(stimsize/Prmt4.s1d)-Prmt4.c2d*erf(stimsize/Prmt4.s2d);
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			L4[i]=1*A4*exp(-(dd*dd*dx2deg2/sigmaL4_2));//rate field of L4
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}
		
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=L4[j]*exp(-(dd*dd*dx2deg2/J.sigmaL4toL23_2));
			}
			H4[i]=dx2deg2*sum;//input current from L4
		}	
	}

	if(Htype==2){
		//INVERSE INPUT is a difference of 2 gaussians in L4 (see Supplement)
		double A4=Prmt4.A1*(erf(stimsize/Prmt4.srho1)-erf(stimsize/Prmt4.srho2));
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			L4[i]=1*A4*(exp(-(dd*dd*dx2deg2/(sigmaL4i_2)))-exp(-(dd*dd*dx2deg2/(sigmaL4bari_2))));//rate field of L4
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=L4[j]*exp(-(dd*dd*dx2deg2/(J.sigmaL4toL23_2)));
			}
			H4[i]=dx2deg2*sum;//input current from L4
		}	
	}
	cout << "rfL4center=" << L4[center] << "\ttotincurrentL4=" << J.L4toPyr*H4[center] << endl;

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

	if(stimsize==0){//this is just to avoid patologies in the Gaussian function
		sigmaLM_2=20;
		sigmaLMi_2=20;
		sigmaLMbari_2=20;
	}

	double *LM = new double[N*N];
	double *LMextra = new double[N*N]; //small excitatory extra input rate field do correct for excitatory input to SOM
	J.LMtoPyr=J.LMtoPyr/(3.141592*J.sigmaLMtoL23_2);


	if(Htype==1){
		//CLASSICAL INPUT is a 2D gaussian field in L4 (see Supplement)
		double ALM=PrmtM.c1d*erf(stimsize/PrmtM.s1d)-PrmtM.c2d*erf(stimsize/PrmtM.s2d);
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			LM[i]=1*ALM*exp(-(dd*dd*dx2deg2/(sigmaLM_2)));//rate field of LM
			//CHANGE CONSTANT HERE FOR CONTRAST OR OPTOGENETIC MODULATIONS
		}

		for(int i=0;i<N*N;i++){
			double sum=0;
			double sum1=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=LM[j]*exp(-(dd*dd*dx2deg2/(J.sigmaLMtoL23_2)));
			}
			HLM[i]=dx2deg2*sum;//input current from LM
		}	
	}
	if(Htype==2){
		//INVERSE INPUT is a difference of 2 gaussians in L4 (see Supplement)
		double ALM1=PrmtM.A1+PrmtM.A2*stimsize;
		double ALM2=ALM1-PrmtM.rho1*(erf(stimsize/PrmtM.srho1)-erf(stimsize/PrmtM.srho2));
		cout<<"ALM1="<<ALM1<<endl;
		cout<<"ALM2="<<ALM2<<endl;
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			LM[i]=1*(ALM1*exp(-(dd*dd*dx2deg2/(sigmaLMi_2)))-ALM2*exp(-(dd*dd*dx2deg2/(sigmaLMbari_2))));//make sure that sigmabari is smaller
			//CHANGE CONSTANT HERE FOR CONTRAST OR OPTOGENETIC MODULATIONS
		}
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=LM[j]*exp(-(dd*dd*dx2deg2/(J.sigmaLMtoL23_2)));
			}
			HLM[i]=dx2deg2*sum;//input current from LM
		}	
	}
	cout << "rfLMcenter=" << LM[center] << "\ttotincurrentLM=" << J.LMtoPyr*HLM[center] << endl;

	delete[] LM;
}



//define the INPUT CURRENT FROM L4 (i.e. a convolution between activity and projections)
void DefineSOMInputs(int Htype, int stimsize){

	//input parameters as defined by the experimental measures (fits/parametrization)
	ReadParamSOM(PrmtS);
	double dd;
	double sigmaS_2=2.*(PrmtS.kdc+PrmtS.kd*stimsize)*(PrmtS.kdc+PrmtS.kd*stimsize);
	double sigmaSi1_2=2.*(PrmtS.kic1+PrmtS.ki*stimsize)*(PrmtS.kic1+PrmtS.ki*stimsize);
	double sigmaSi2_2=2.*(PrmtS.kic2+PrmtS.ki*stimsize)*(PrmtS.kic2+PrmtS.ki*stimsize);
	J.SOMtoPyr=J.SOMtoPyr/(3.141592*J.sigmaSOMtoPyr_2);//Gaussian connectivity is normalized

	//define spatial dependence of the FeedForward current from layer 4
	double *SOM = new double[N*N];

	if(Htype==1){
		//CLASSICAL INPUT is a 2D gaussian field in L4 (see Supplement)
		double AS=PrmtS.c1d*erf(stimsize/PrmtS.s1d)-PrmtS.c2d*erf(stimsize/PrmtS.s2d);
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			SOM[i]=1*AS*exp(-(dd*dd*dx2deg2/sigmaS_2));//rate field of SOM
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}
		
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=SOM[j]*exp(-(dd*dd*dx2deg2/J.sigmaSOMtoPyr_2));
			}
			HS[i]=dx2deg2*sum;//input current from SOM
		}	
	}

	if(Htype==2){
		//INVERSE INPUT is a difference of 2 gaussians in L4 (see Supplement)
		//rRECI[x_, y_, s_] :=  (r0S[s] Exp[-(x^2 + y^2)/(2 bar\[Sigma]RECI1[s])] - \[Rho]S[s] Exp[-(x^2 + y^2)/(2 bar\[Sigma]RECI2[s])]);
		//r0S[s_] := A1S (Erf[s/sA1S] - A2S Erf[s/sA2S]);
		//\[Rho]S[s_] := \[Rho]1S (Erf[(s - 10)/s\[Rho]1S] - \[Rho]2S Erf[(s - 10)/s\[Rho]2S]);
		double r0S=PrmtS.A1*(erf(stimsize/PrmtS.sa1)- PrmtS.A2*erf(stimsize/PrmtS.sa2));//check inverse
		double rhoS=PrmtS.rho1*(erf((stimsize-10)/PrmtS.srho1)- PrmtS.rho2*erf((stimsize-10)/PrmtS.srho2));//check inverse
		for(int i=0;i<N*N;i++){
			dd=finddistance(i,center);
			SOM[i]=1*r0S*exp(-(dd*dd*dx2deg2/(sigmaSi1_2)))-rhoS*exp(-(dd*dd*dx2deg2/(sigmaSi2_2)));//rate field of SOM
			//CHANGE CONSTANT HERE FOR CONTRAST MODULATIONS
		}
		for(int i=0;i<N*N;i++){
			double sum=0;
			for (int j=0;j<N*N;j++){
				dd=finddistance(i,j);
				sum+=SOM[j]*exp(-(dd*dd*dx2deg2/(J.sigmaSOMtoPyr_2)));
			}
			HS[i]=dx2deg2*sum;//input current from L4
		}	
	}
	cout << "rfSOMcenter=" << SOM[center] << "\ttotincurrentSOM=" << J.SOMtoPyr*HS[center] << endl;

	delete[] SOM;
}



void Integrate2DPowerLawIO(double dt){

	int i,j;
	double TotInputPyr,RecInputPyr,RPyrss;

	//Euler rule//
	//update L23
	for(i=0;i<N*N;i++){
		RecInputPyr=0;

		for(j=0;j<N*N;j++)
			RecInputPyr+=WPyrPyr[i][j]*FiringRatesL3.Pyr[j];

		TotInputPyr=J.L4toPyr*H4[i]+dx2deg2*RecInputPyr+J.LMtoPyr*HLM[i]+J.SOMtoPyr*HS[i];

		RPyrss=max(0.,TotInputPyr);
		RPyrss=k*pow(RPyrss,n);
		
		FiringRatesL3.Pyr[i]+=dt*(1/tauPyr)*(-FiringRatesL3.Pyr[i]+RPyrss);
	}
}


double finddistance(int p1,int p2){
	//function to find distance of two points in a 2D lattice of linear size N

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
	//cout<<"p1= "<<p1<<";xp1= "<<xp1<<";yp1= "<<yp1<<endl<<"p2= "<<p2<<";xp2= "<<xp2<<";yp2= "<<yp2<<endl<<"d= "<<d<<endl;
	return d;
}


void InitializeConnections(jij &J){


	//initializes matrices WXY[i][j] i,j=1..N such that interactions depend on the distance(i,j).
	//a unit interacts with neighbors (similar retinotopic locations).
	int i,j;
	double d, d2;

	//here Im normalizing Jij!
	int center=round(N*(N-1)/2);
	
	for(i=0;i<N*N;i++){
		for(j=0;j<N*N;j++){//i want WEE[k][k]=JEE;
			d=finddistance(i,j);
			d2=(double)d*d*dx2deg2;
			WPyrPyr[i][j]=(J.PyrPyr*exp(-d2/J.sigma2PyrPyr2))/(3.141592*J.sigma2PyrPyr2);
		}
	}
	//cout << WPyrPyr[0][0] << endl << endl;

}



void CompareWithAnalytics(void){

	double RecInputPyr;

	RecInputPyr=0;

	for(int j=0;j<N*N;j++)
		RecInputPyr+=WPyrPyr[center][j]*FiringRatesL3.Pyr[j];

	cout << "Pyrinput=" << dx2deg2*RecInputPyr << "\tL4input" << J.L4toPyr*H4[center] << "\tLMinput" << J.LMtoPyr*HLM[center] << "\tSOMinput" << J.SOMtoPyr*HS[center] << endl;

}




