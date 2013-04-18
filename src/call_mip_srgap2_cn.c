/*
 * Xander Nuttle
 *
 * call_mip_srgap2_cn.c
 *
 * Generates paralog-specific SRGAP2 copy number calls from MIP paralog-specific
 * read count data, using a maximum likelihood approach and dynamic programming.
 *
 * Call: ./call_mip_srgap2_cn SRGAP2_miptargets_file SRGAP2_mip_read_counts_file output_file_base_name
 * Example call: ./call_mip_srgap2_cn SRGAP2_final.miptargets pos_ctrl_expt_SRGAP2_final.mipcounts pos_ctrl_expt
 *
 * NOTE: The mipcounts file must contain data for only MIPs in the miptargets
 * file!
 *
 * This program takes in (1) a file containing information on SRGAP2 MIP targets
 * and (2) a file containing paralog-specific SRGAP2 read counts, outputting (1)
 * a file containing SRGAP2 copy number genotypes for each individual included
 * in the paralog-specific SRGAP2 read counts file. Input file (2) should
 * contain data corresponding only to MIPs in input file (1); that is, the data
 * in file 2 should be generated by running mip_analysis.c with file 1 as an
 * input.  Additionally, all MIPs in file (1) must be sorted according to master
 * sequence target coordinate, smallest to largest.
 *
 * This program generates genotype calls using a maximum likelihood approach
 * with dynamic programming. For each individual, likelihoods of observing the
 * count data for each MIP are calculated under 400 different possible
 * underlying SRGAP2 paralog-specific copy number states, where SRGAP2A and
 * SRGAP2C can have copy number 0-3 and SRGAP2B and SRGAP2D can have copy number
 * 0-4 (400 combinations = 4*5*4*5). Then, for each individual, these
 * likelihoods are used to construct a weighted directed acyclic graph (WDAG),
 * with prior probabilities based on observed genotype data from previous
 * experiments incorporated into the likelihoods for the first MIP count
 * data. As the graph is constructed, highest scoring paths allowing 1 and 2
 * transitions between copy number states are tracked. Transitions are
 * restricted to duplications or deletions affecting 1 paralog only or gene
 * conversion events, where the copy numbers of two paralogs change, but the
 * total number of SRGAP2 copies remains constant. Placing these restrictions on
 * transitions reflects the fact that true biological events should fall into
 * one of these categories (single-paralog-affecting duplicaion/deletion or gene
 * conversion). The program ultimately identifies the highest-scoring (most
 * likely) paths through the likelihood graph allowing 0, 1, and 2 transitions
 * and their corresponding scores (likelihoods). A heuristic is used to assess
 * increases in likelihood of the 1-transition and 2-transition paths compared
 * to the 0 transition path and determine whether they signal an internal event
 * (duplication, deletion, or conversion) and warrant calling the genotype for
 * an individual as complex. In most cases, the scores of the 1-transition and
 * 2-transition paths will not be substantially higher than that of the
 * 0-transition path, and the genotype for an individual will be called as
 * simple (a single copy number state across the entirety of the gene).
 */

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

//MAGIC NUMBERS
//static const int NUM_PLOGS=4; //number of paralogs in gene family
//static const int SRGAP2A_MAX_CN=3; //maximum copy number genotype for SRGAP2A
//static const int SRGAP2B_MAX_CN=4; //maximum copy number genotype for SRGAP2B
//static const int SRGAP2C_MAX_CN=3; //maximum copy number genotype for SRGAP2C
//static const int SRGAP2D_MAX_CN=4; //maximum copy number genotype for SRGAP2D
#define NUM_PLOGS 4
#define SRGAP2A_MAX_CN 3
#define SRGAP2B_MAX_CN 4
#define SRGAP2C_MAX_CN 3
#define SRGAP2D_MAX_CN 4
//static const int MIN_LIKELIHOOD=-30; //parameter setting the minimum likelihood value for a single data point
//static const int NUM_CN_STATES=(SRGAP2A_MAX_CN+1)*(SRGAP2B_MAX_CN+1)*(SRGAP2C_MAX_CN+1)*(SRGAP2D_MAX_CN+1); //number of distinct possible copy number states
#define MIN_LIKELIHOOD -30
#define NUM_CN_STATES ((SRGAP2A_MAX_CN+1)*(SRGAP2B_MAX_CN+1)*(SRGAP2C_MAX_CN+1)*(SRGAP2D_MAX_CN+1))
#define SRGAP2A_FREQ_0 exp(MIN_LIKELIHOOD)
#define SRGAP2C_FREQ_0 exp(MIN_LIKELIHOOD)
//static const double SRGAP2A_FREQ_0=exp(MIN_LIKELIHOOD); //frequency of a SRGAP2A copy number genotype of 0
static const double SRGAP2A_FREQ_1=(double)3/(double)28153; //frequency of a SRGAP2A copy number genotype of 1
static const double SRGAP2A_FREQ_2=(double)28147/(double)28153; //frequency of a SRGAP2A copy number genotype of 2
static const double SRGAP2A_FREQ_3=(double)3/(double)28153; //frequency of a SRGAP2A copy number genotype of 3
static const double SRGAP2B_FREQ_0=(double)3/(double)661; //frequency of a SRGAP2B copy number genotype of 0
static const double SRGAP2B_FREQ_1=(double)33/(double)661; //frequency of a SRGAP2B copy number genotype of 1
static const double SRGAP2B_FREQ_2=(double)575/(double)661; //frequency of a SRGAP2B copy number genotype of 2
static const double SRGAP2B_FREQ_3=(double)49/(double)661; //frequency of a SRGAP2B copy number genotype of 3
static const double SRGAP2B_FREQ_4=(double)1/(double)661; //frequency of a SRGAP2B copy number genotype of 4
//static const double SRGAP2C_FREQ_0=exp(MIN_LIKELIHOOD); //frequency of a SRGAP2C copy number genotype of 0
static const double SRGAP2C_FREQ_1=(double)2/(double)7140; //frequency of a SRGAP2C copy number genotype of 1
static const double SRGAP2C_FREQ_2=(double)7134/(double)7140; //frequency of a SRGAP2C copy number genotype of 2
static const double SRGAP2C_FREQ_3=(double)4/(double)7140; //frequency of a SRGAP2C copy number genotype of 3
static const double SRGAP2D_FREQ_0=(double)3/(double)47; //frequency of a SRGAP2D copy number genotype of 0
static const double SRGAP2D_FREQ_1=(double)9/(double)47; //frequency of a SRGAP2D copy number genotype of 1
static const double SRGAP2D_FREQ_2=(double)28/(double)47; //frequency of a SRGAP2D copy number genotype of 2
static const double SRGAP2D_FREQ_3=(double)6/(double)47; //frequency of a SRGAP2D copy number genotype of 3
static const double SRGAP2D_FREQ_4=(double)1/(double)47; //frequency of a SRGAP2D copy number genotype of 4
static const long SRGAP2D_DEL_START=105947; //master sequence base 1 coordinate of SRGAP2D deletion region start
static const long SRGAP2D_DEL_END=213356; //master sequence base 1 coordinate of SRGAP2D deletion region end
static const int MIN_LIKELIHOOD_DIFF=40; //see below for detailed description of this heuristic
static const int MIN_MIPS_IN_CN_STATE=5; //see below for detailed description of this heuristic

int main(int argc,char*argv[])
{
	//get information about number of mip targets designed for copy number genotyping
	FILE*miptargetsfile;
	fpos_t pos;
	char miptype;
	char dummy[501];
	long num_mip_targets=0,num_exon_mips=0;
	miptargetsfile=fopen(*(argv+1),"r");
	fgetpos(miptargetsfile,&pos);
	while(fscanf(miptargetsfile,"%s %s %s %c %s %s %s %s %s",dummy,dummy,dummy,&miptype,dummy,dummy,dummy,dummy,dummy)==9)
	{
		if(miptype!='E')
		{
			num_mip_targets++;
		}
		else
		{
			num_exon_mips++;
		}
	}
	fsetpos(miptargetsfile,&pos);

	//get and store information about SRGAP2 master sequence target coordinate and specificity for each MIP
	char spec_string[NUM_PLOGS];
	char spec[2]; // a string to hold either 0 or 1 (a MIP either has or lacks specificity for a given paralog)
	long start,end;
	long target_coords[num_mip_targets];
	long specificities[num_mip_targets][NUM_PLOGS];
	long i=0,j;
	spec[1]='\0';
	while(fscanf(miptargetsfile,"%s %ld %c %ld %s %c %s %s %s %s %s",dummy,&start,&miptype,&end,dummy,&miptype,dummy,dummy,spec_string,dummy,dummy)==11)
  {
		if(miptype!='E')
    {
			target_coords[i]=(start+end)/2;
			for(j=0;j<NUM_PLOGS;j++)
			{
				spec[0]=spec_string[j];
				specificities[i][j]=strtol(spec,NULL,10);
			}
    	i++;
		}
  }
	fclose(miptargetsfile);

	//generate vector of possible SRGAP2 copy number states and vector of ln(prior probabilities) corresponding to each copy number state
	//assume independence of paralog-specific copy number genotypes
	int A,B,C,D;
  int copy_states[NUM_CN_STATES][NUM_PLOGS];
  double priors[NUM_CN_STATES]; //vector of prior probabilities for each copy number state
		//priors_N gives the probabilities of observing a copy number state of 0,1,2,3,... for SRGAP2N
		//these probabilities are based on empirical data from all SRGAP2 duplication/deletion screening and genotyping efforts
		//extremely unlikely states (never observed so far but hypothetically possible) are arbitrarily assigned a value with a resulting natural logarithm of -30
  double priors_A[SRGAP2A_MAX_CN+1]={SRGAP2A_FREQ_0,SRGAP2A_FREQ_1,SRGAP2A_FREQ_2,SRGAP2A_FREQ_3};
  double priors_B[SRGAP2B_MAX_CN+1]={SRGAP2B_FREQ_0,SRGAP2B_FREQ_1,SRGAP2B_FREQ_2,SRGAP2B_FREQ_3,SRGAP2B_FREQ_4};
  double priors_C[SRGAP2C_MAX_CN+1]={SRGAP2C_FREQ_0,SRGAP2C_FREQ_1,SRGAP2C_FREQ_2,SRGAP2C_FREQ_3}; //the 3 Signature Genomics cases genotyped by MIPs are included in denominator
  double priors_D[SRGAP2D_MAX_CN+1]={SRGAP2D_FREQ_0,SRGAP2D_FREQ_1,SRGAP2D_FREQ_2,SRGAP2D_FREQ_3,SRGAP2D_FREQ_4}; //all SRGAP2D genotyping is from manual inspection of MIP data, Troina individual omitted
  i=0;
	for(A=0;A<=SRGAP2A_MAX_CN;A++)
  {
    for(B=0;B<=SRGAP2B_MAX_CN;B++)
    {
      for(C=0;C<=SRGAP2C_MAX_CN;C++)
      {
        for(D=0;D<=SRGAP2D_MAX_CN;D++)
        {
					copy_states[i][0]=A;
          copy_states[i][1]=B;
          copy_states[i][2]=C;
          copy_states[i][3]=D;
          priors[i]=log(priors_A[A])+log(priors_B[B])+log(priors_C[C])+log(priors_D[D]);
          i++;
        }
      }
    }
  }

	//set up node structure for dynamic programming to find maximum likelihood path through graph
	//transitions to allow = 2*(number of internal events to allow), because each internal event can maximally contribute two edges in copy number state corresponding to 2 transitions
	//right now, the program allows for 1 internal event (testing to see if this will still pick up > 1 internal event, in which case allowing for more internal events will not be necessary)
	struct node
	{
		long state;
		double max0;
		double max1;
		double max2;
		struct node*prev_node_1trans;
		struct node*prev_node_2trans;
		long index;
	};

	//allocate storage for likelihood graph
	struct node*likelihood_graph;
	likelihood_graph=(struct node*)malloc(num_mip_targets*NUM_CN_STATES*sizeof(struct node));

	//initialize output files
	FILE*out,*out2,*out3;
  char output_base_name[51];
  char output_file_extension[9]=".cncalls";
  char ext2[12]=".compevents";
	char ext3[13]=".simplecalls";
	strncpy(output_base_name,*(argv+3),42);
  strcat(output_base_name,output_file_extension);
  out=fopen(output_base_name,"w");
	for(i=0;i<51;i++)
	{
		output_base_name[i]='\0';
	}
	strncpy(output_base_name,*(argv+3),39);
  strcat(output_base_name,ext2);
  out2=fopen(output_base_name,"w");
	for(i=0;i<51;i++)
  {
    output_base_name[i]='\0';
  }
  strncpy(output_base_name,*(argv+3),38);
  strcat(output_base_name,ext3);
  out3=fopen(output_base_name,"w");
	fprintf(out3,"Individual\tSRGAP2A_CN\tSRGAP2B_CN\tSRGAP2C_CN\tSRGAP2D_CN\tPossible_Complex_CN_Genotype\n"); //THIS ISN'T A MAGIC NUMBER BUT IS A CONSTANT OF SORTS

	//for each individual, read in counts, calculate and store individual likelihoods of data for each MIP under each copy number state
	FILE*countsfile;
	char individual[31];
	long counts[NUM_PLOGS+1];
	long indiv_counts[num_mip_targets][NUM_PLOGS+1];
	long coord,mip_coord;
	double probs[NUM_PLOGS+1];
	double L;
	long k;
	double mip_likelihoods[num_mip_targets][NUM_CN_STATES];
	individual[30]='\0';
	countsfile=fopen(*(argv+2),"r");
	while(getc(countsfile)!='\n')
    continue;
	fgetpos(countsfile,&pos);
	while(fscanf(countsfile,"%s %s %ld %c %ld %ld %ld %ld %ld",individual,dummy,&coord,&miptype,&(counts[0]),&(counts[1]),&(counts[2]),&(counts[3]),&(counts[4]))==9) //while there remains data to process
	{
		fsetpos(countsfile,&pos);
		j=0;
		for(i=0;i<(num_mip_targets+num_exon_mips);i++)
		{
			fscanf(countsfile,"%s %s %ld %c %ld %ld %ld %ld %ld",individual,dummy,&mip_coord,&miptype,&(counts[0]),&(counts[1]),&(counts[2]),&(counts[3]),&(counts[4]));
			if(miptype!='E')
			{
				for(k=0;k<(NUM_PLOGS+1);k++)
				{
					indiv_counts[j][k]=counts[k];
				}
				j++;
			}
		}
		//calculate likelihoods of observed counts at each MIP under each possible copy number state
		for(i=0;i<num_mip_targets;i++)
		{
			for(j=0;j<NUM_CN_STATES;j++)
			{
				probs[NUM_PLOGS]=0.0;
				if((target_coords[i]<SRGAP2D_DEL_START)||(target_coords[i]>SRGAP2D_DEL_END)) //MIP not in SRGAP2D deletion region
				{
					for(k=0;k<NUM_PLOGS;k++)
					{
						if(specificities[i][k]==1)
						{
							probs[k]=(double)(copy_states[j][k])/(double)(copy_states[j][0]+copy_states[j][1]+copy_states[j][2]+copy_states[j][3]);
						}
						else
						{
							probs[k]=0.0;
							probs[NUM_PLOGS]+=(double)(copy_states[j][k])/(double)(copy_states[j][0]+copy_states[j][1]+copy_states[j][2]+copy_states[j][3]);
						}
					}
				}
				else
				{
					probs[NUM_PLOGS-1]=0.0;
					for(k=0;k<(NUM_PLOGS-1);k++)
          {
            if(specificities[i][k]==1)
            {
              probs[k]=(double)(copy_states[j][k])/(double)(copy_states[j][0]+copy_states[j][1]+copy_states[j][2]);
            }
            else
            {
              probs[k]=0.0;
              probs[NUM_PLOGS]+=(double)(copy_states[j][k])/(double)(copy_states[j][0]+copy_states[j][1]+copy_states[j][2]);
            }
          }
				}
				const unsigned int countvec[NUM_PLOGS+1]={indiv_counts[i][0],indiv_counts[i][1],indiv_counts[i][2],indiv_counts[i][3],indiv_counts[i][4]};
				L=gsl_ran_multinomial_lnpdf(NUM_PLOGS+1,probs,countvec);
				if(L<MIN_LIKELIHOOD)
				{
					L=MIN_LIKELIHOOD; //minimum log likelihood value for 1 MIP probe arbitrarily assigned to -30
				}
				mip_likelihoods[i][j]=L; //[mip_target][copy_state]
			}
		}

		//use dynamic programming to determine whether the data support more than 1 copy number state across SRGAP2
		double max_max0=-DBL_MAX,max_max1=-DBL_MAX,max_max2=-DBL_MAX;
    long maxindex_max0,maxindex_max1,maxindex_max2;
			//setup first column and calculate maximum likelihood value and index of corresponding node
		for(j=0;j<NUM_CN_STATES;j++)
    {
      likelihood_graph[j].index=j;
			likelihood_graph[j].state=j;
			likelihood_graph[j].max0=mip_likelihoods[0][j]+priors[j];
      likelihood_graph[j].max1=-DBL_MAX;
      likelihood_graph[j].max2=-DBL_MAX;
			likelihood_graph[j].prev_node_1trans=NULL;
			likelihood_graph[j].prev_node_2trans=NULL;
    }
			//setup all other verticies
		for(j=NUM_CN_STATES;j<(NUM_CN_STATES*num_mip_targets);j++)
    {
      if((j%NUM_CN_STATES)==0)
      {
        max_max0=-DBL_MAX;
        max_max1=-DBL_MAX;
        for((k=j-NUM_CN_STATES);k<j;k++)
        {
          if(likelihood_graph[k].max0>max_max0)
          {
            max_max0=likelihood_graph[k].max0;
						maxindex_max0=k;
          }
          if(likelihood_graph[k].max1>max_max1)
          {
            max_max1=likelihood_graph[k].max1;
						maxindex_max1=k;
          }
        }
      }
      likelihood_graph[j].index=j;
			likelihood_graph[j].state=j%NUM_CN_STATES;
			likelihood_graph[j].max0=mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES]+likelihood_graph[j-NUM_CN_STATES].max0;
      if(max_max0>likelihood_graph[j-NUM_CN_STATES].max1)
      {
        likelihood_graph[j].max1=max_max0+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
				likelihood_graph[j].prev_node_1trans=&(likelihood_graph[maxindex_max0]);
      }
      else
      {
        likelihood_graph[j].max1=likelihood_graph[j-NUM_CN_STATES].max1+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
      	likelihood_graph[j].prev_node_1trans=&(likelihood_graph[j-NUM_CN_STATES]);
			}
      if(max_max1>likelihood_graph[j-NUM_CN_STATES].max2)
      {
        likelihood_graph[j].max2=max_max1+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
				likelihood_graph[j].prev_node_2trans=&(likelihood_graph[maxindex_max1]);
      }
      else
      {
        likelihood_graph[j].max2=likelihood_graph[j-NUM_CN_STATES].max2+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
				likelihood_graph[j].prev_node_2trans=&(likelihood_graph[j-NUM_CN_STATES]);
      }
		}
			//determine maximum values of final likelihoods allowing 0, 1, and 2 transitions and use to determine whether there may be multiple underlying copy number states across SRGAP2
		max_max0=-DBL_MAX;
    max_max1=-DBL_MAX;
    max_max2=-DBL_MAX;
    for(j=(NUM_CN_STATES*num_mip_targets-NUM_CN_STATES);j<(NUM_CN_STATES*num_mip_targets);j++)
    {
      if(likelihood_graph[j].max0>max_max0)
      {
        max_max0=likelihood_graph[j].max0;
  			maxindex_max0=j%NUM_CN_STATES; //maximally likely copy number state assuming no internal events
			}
			if(likelihood_graph[j].max1>max_max1)
      {
        max_max1=likelihood_graph[j].max1;
				maxindex_max1=j;
      }
      if(likelihood_graph[j].max2>max_max2)
      {
        max_max2=likelihood_graph[j].max2;
      	maxindex_max2=j;
			}
		}
		fprintf(out,"Individual: %s\nL0=%lf L1-L0=%lf L2-L0=%lf\n",individual,max_max0,max_max1-max_max0,max_max2-max_max0);
		fprintf(out,"Genotype (assuming no internal events): %d%d%d%d\n",copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],copy_states[maxindex_max0][2],copy_states[maxindex_max0][3]);
		int evidence_for_1T=0,evidence_for_2T=0;
		long num_mips_state1=0,num_mips_state2=0,num_mips_state3=0;
    struct node*current;
    long oldstate;
    int trans_made=0;
			//
			//Heuristics for predicitng 1 copy number state transition
			//	Likelihood_1_transition - Likelihood_0_transitions > 40 (MIN_LIKELIHOOD_DIFF)
			//	At least 5 MIPs in each copy number state (MIN_MIPS_IN_CN_STATE)
			//	Likelihood_2_transitions - Likelihood_1_transition <= 40
			//
			//Heuristics for prediciting 2 copy number state transitions
			//	Likelihood_2_transitions - Likelihood_0_transitions > 40
			//	At least 5 MIPs in each copy number state
			//	Likelihood_2_transitions - Likelihood_1_transition > 40
			//
		if((max_max1-max_max0)>MIN_LIKELIHOOD_DIFF)
		{
			//trace back to determine number of MIPs in each state
			//1 transition
			current=&(likelihood_graph[maxindex_max1]);
    	while(1)
    	{
      	num_mips_state2++;
				oldstate=current->state;
      	current=current->prev_node_1trans;
      	if((current->state)!=oldstate)
      	{
        	break;
      	}
    	}
			num_mips_state1=num_mip_targets-num_mips_state2;
			if((num_mips_state1>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state2>=MIN_MIPS_IN_CN_STATE))
			{
				evidence_for_1T=1;
			}
		}
		if((max_max2-max_max0)>MIN_LIKELIHOOD_DIFF)
		{
			//2 transitions
			num_mips_state1=0;
			num_mips_state2=0;
			current=&(likelihood_graph[maxindex_max2]);
			while(trans_made<2)
			{
				oldstate=current->state;
				if(trans_made==0)
				{
					num_mips_state3++;
					current=current->prev_node_2trans;
				}
				else
				{
					num_mips_state2++;
					current=current->prev_node_1trans;
				}
				if((current->state)!=oldstate)
				{
					trans_made++;
				}
			}
			num_mips_state1=num_mip_targets-num_mips_state2-num_mips_state3;
			if((num_mips_state1>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state2>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state3>=MIN_MIPS_IN_CN_STATE))
			{
				if(evidence_for_1T)
				{
					if(max_max2-max_max1>MIN_LIKELIHOOD_DIFF)
					{
						evidence_for_2T=1;
						evidence_for_1T=0;
					}
				}
				else
				{
					evidence_for_2T=1;
					evidence_for_1T=0;
				}
			}
		}


		//if there is evidence for possible multiple copy number states across SRGAP2, use dynamic programming to find maximum likelihood paths through likelihood graph, allowing up to 2 transitions (1 internal event)
		int num_paralog_cn_changes,total_copies_prev,total_copies;
		if(evidence_for_1T||evidence_for_2T)
		{
			//alert user
			printf("\n%s may have a complex genotype!!! Graphically view the data for %s to ensure genotyping accuracy.\n",individual,individual);
			//print to 2nd output file so you can easily visualize only these individuals
			fprintf(out2,"%s\t",individual);
			fprintf(out2,"L0=%lf L1-L0=%lf L2-L0=%lf\t",max_max0,max_max1-max_max0,max_max2-max_max0);
    	fprintf(out2,"Genotype (assuming no internal events): %d%d%d%d\t",copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],copy_states[maxindex_max0][2],copy_states[maxindex_max0][3]);
			fprintf(out3,"%s\t%d\t%d\t%d\t%d\tYES\n",individual,copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],copy_states[maxindex_max0][2],copy_states[maxindex_max0][3]);
			//setup first column and calculate maximum likelihood value and index of corresponding node
			for(j=0;j<NUM_CN_STATES;j++)
			{
				likelihood_graph[j].state=j;
				likelihood_graph[j].index=j;
				likelihood_graph[j].max0=mip_likelihoods[0][j]+priors[j];
				likelihood_graph[j].max1=-DBL_MAX;
				likelihood_graph[j].max2=-DBL_MAX;
				likelihood_graph[j].prev_node_1trans=NULL;
				likelihood_graph[j].prev_node_2trans=NULL;
			}
			//setup all other verticies
			for(j=NUM_CN_STATES;j<(NUM_CN_STATES*num_mip_targets);j++)
			{
      	max_max0=-DBL_MAX;
      	max_max1=-DBL_MAX;
				maxindex_max0=0;
				maxindex_max1=0;
				for((k=j-NUM_CN_STATES);k<j;k++)
      	{
        	num_paralog_cn_changes=0;
					total_copies_prev=0;
					total_copies=0;
					for(i=0;i<4;i++)
					{
						if(copy_states[k%NUM_CN_STATES][i]!=copy_states[j%NUM_CN_STATES][i])
						{
							num_paralog_cn_changes++;
							total_copies_prev+=copy_states[k%NUM_CN_STATES][i];
							total_copies+=copy_states[j%NUM_CN_STATES][i];
						}
					}
					if((num_paralog_cn_changes==2)&&(total_copies_prev==total_copies)) //signature of gene conversion, still a single event
        	{
          	num_paralog_cn_changes=1;
        	}
					if((likelihood_graph[k].max0>max_max0)&&(num_paralog_cn_changes==1)) //only allow a single event
        	{
          	max_max0=likelihood_graph[k].max0;
          	maxindex_max0=k;
        	}
					if((likelihood_graph[k].max1>max_max1)&&(num_paralog_cn_changes==1)) //only allow a single event
					{
						max_max1=likelihood_graph[k].max1;
          	maxindex_max1=k;
					}
      	}
				likelihood_graph[j].state=j%NUM_CN_STATES;
				likelihood_graph[j].index=j;
				likelihood_graph[j].max0=mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES]+likelihood_graph[j-NUM_CN_STATES].max0;
				if(max_max0>likelihood_graph[j-NUM_CN_STATES].max1)
				{
					likelihood_graph[j].max1=max_max0+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
					likelihood_graph[j].prev_node_1trans=&(likelihood_graph[maxindex_max0]);
				}
				else
				{
					likelihood_graph[j].max1=likelihood_graph[j-NUM_CN_STATES].max1+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
					likelihood_graph[j].prev_node_1trans=&(likelihood_graph[j-NUM_CN_STATES]);
				}
				if(max_max1>likelihood_graph[j-NUM_CN_STATES].max2)
				{
					likelihood_graph[j].max2=max_max1+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
        	likelihood_graph[j].prev_node_2trans=&(likelihood_graph[maxindex_max1]);
				}
				else
				{
					likelihood_graph[j].max2=likelihood_graph[j-NUM_CN_STATES].max2+mip_likelihoods[j/NUM_CN_STATES][j%NUM_CN_STATES];
        	likelihood_graph[j].prev_node_2trans=&(likelihood_graph[j-NUM_CN_STATES]);
				}
			}
			//determine maximum values of final likelihoods allowing 0, 1, and 2 transitions
			max_max0=-DBL_MAX;
			max_max1=-DBL_MAX;
			max_max2=-DBL_MAX;
			maxindex_max2=0;
			for(j=(NUM_CN_STATES*num_mip_targets-NUM_CN_STATES);j<(NUM_CN_STATES*num_mip_targets);j++)
			{
				if(likelihood_graph[j].max0>max_max0)
				{
					max_max0=likelihood_graph[j].max0;
					maxindex_max0=j%NUM_CN_STATES; //maximally likely copy number state assuming no internal events
				}
				if(likelihood_graph[j].max1>max_max1)
				{
					max_max1=likelihood_graph[j].max1;
					maxindex_max1=j;
				}
				if(likelihood_graph[j].max2>max_max2)
				{
					max_max2=likelihood_graph[j].max2;
        	maxindex_max2=j;
				}
			}
			//trace back to determine copy number states at each MIP for 1 and 2 transition scenarios
			struct node*current;
				//1 transition
			long first_mip_state2_1trans,state1_1trans,oldstate,oldindex;
			current=&(likelihood_graph[maxindex_max1]);
			while((current->index)>=NUM_CN_STATES)
			{
				oldstate=current->state;
				oldindex=current->index;
				current=current->prev_node_1trans;
				if((current->state)!=oldstate)
				{
					first_mip_state2_1trans=target_coords[oldindex/NUM_CN_STATES];
					state1_1trans=current->state;
					break;
				}
			}
				//2 transitions
			long first_mip_state2_2trans,last_mip_state2_2trans,state2_2trans,state1_2trans;
			current=&(likelihood_graph[maxindex_max2]);
    	int trans_made=0;
			while(trans_made<2)
			{
				oldstate=current->state;
				oldindex=current->index;
				if(trans_made==0)
				{
					current=current->prev_node_2trans;
				}
				else
				{
					current=current->prev_node_1trans;
				}
				if((current->state)!=oldstate)
				{
					if(trans_made==0)
					{
						state2_2trans=current->state;
						last_mip_state2_2trans=target_coords[(current->index)/NUM_CN_STATES];
						trans_made++;
					}
					else
					{
						state1_2trans=current->state;
						first_mip_state2_2trans=target_coords[oldindex/NUM_CN_STATES];
						trans_made++;
					}
				}
			}
			if(evidence_for_1T)
			{
				fprintf(out,"Genotype (assuming 1 transition): %d%d%d%d %d%d%d%d, First MIP in State 2: %ld\n\n",copy_states[state1_1trans][0],copy_states[state1_1trans][1],copy_states[state1_1trans][2],copy_states[state1_1trans][3],copy_states[maxindex_max1%NUM_CN_STATES][0],copy_states[maxindex_max1%NUM_CN_STATES][1],copy_states[maxindex_max1%NUM_CN_STATES][2],copy_states[maxindex_max1%NUM_CN_STATES][3],first_mip_state2_1trans);
			fprintf(out2,"Genotype (assuming 1 transition): %d%d%d%d %d%d%d%d, First MIP in State 2: %ld\n",copy_states[state1_1trans][0],copy_states[state1_1trans][1],copy_states[state1_1trans][2],copy_states[state1_1trans][3],copy_states[maxindex_max1%NUM_CN_STATES][0],copy_states[maxindex_max1%NUM_CN_STATES][1],copy_states[maxindex_max1%NUM_CN_STATES][2],copy_states[maxindex_max1%NUM_CN_STATES][3],first_mip_state2_1trans);
      }
			else
			{
				fprintf(out,"Genotype (assuming 2 transitions): %d%d%d%d %d%d%d%d %d%d%d%d, First MIP in State 2: %ld, Last MIP in State 2: %ld\n\n",copy_states[state1_2trans][0],copy_states[state1_2trans][1],copy_states[state1_2trans][2],copy_states[state1_2trans][3],copy_states[state2_2trans][0],copy_states[state2_2trans][1],copy_states[state2_2trans][2],copy_states[state2_2trans][3],copy_states[maxindex_max2%NUM_CN_STATES][0],copy_states[maxindex_max2%NUM_CN_STATES][1],copy_states[maxindex_max2%NUM_CN_STATES][2],copy_states[maxindex_max2%NUM_CN_STATES][3],first_mip_state2_2trans,last_mip_state2_2trans);
			fprintf(out2,"Genotype (assuming 2 transitions): %d%d%d%d %d%d%d%d %d%d%d%d, First MIP in State 2: %ld, Last MIP in State 2: %ld\n",copy_states[state1_2trans][0],copy_states[state1_2trans][1],copy_states[state1_2trans][2],copy_states[state1_2trans][3],copy_states[state2_2trans][0],copy_states[state2_2trans][1],copy_states[state2_2trans][2],copy_states[state2_2trans][3],copy_states[maxindex_max2%NUM_CN_STATES][0],copy_states[maxindex_max2%NUM_CN_STATES][1],copy_states[maxindex_max2%NUM_CN_STATES][2],copy_states[maxindex_max2%NUM_CN_STATES][3],first_mip_state2_2trans,last_mip_state2_2trans);
      }
		}
		else
		{
			fprintf(out,"\n");
			fprintf(out3,"%s\t%d\t%d\t%d\t%d\tNO\n",individual,copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],copy_states[maxindex_max0][2],copy_states[maxindex_max0][3]);
			if((copy_states[maxindex_max0][0]!=2)||copy_states[maxindex_max0][2]!=2)
			{
				printf("\nSRGAP2A/C EVENT DETECTED IN %s!!! Genotype: %d%d%d%d\n",individual,copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],copy_states[maxindex_max0][2],copy_states[maxindex_max0][3]);
			}
		}

		//prepare for next individual
		fgetpos(countsfile,&pos);
	}

	//clean up and exit
	free(likelihood_graph);
	fclose(out);
	fclose(out2);
	fclose(out3);
	return 0;
}

