//Xander Nuttle
//call_mip_rh_cn.c
//Generates paralog-specific RH copy number calls from MIP paralog-specific read count data, using a maximum likelihood approach and dynamic programming.
//Call: ./call_mip_rh_cn RH_miptargets_file RH_mip_read_counts_file output_file_base_name
//Example call: ./call_mip_rh_cn RH_gene_only.miptargets pos_ctrl_expt_RH.mipcounts pos_ctrl_expt_RH
//
//THE MIPCOUNTS FILE MUST CONTAIN DATA FOR ONLY MIPS IN THE MIPTARGETS FILE!!!
//
//This program takes in (1) a file containing information on RH MIP targets and (2) a file containing paralog-specific RH read counts, outputting (1)
//a file containing RH copy number genotypes and LOD score confidence values for each individual included in the paralog-specific RH read counts file.
//Input file (2) should contain data corresponding only to MIPs in input file (1); that is, the data in file 2 should be generated by running the
//paralog-specific read counting program with file 1 as an input. Additionally, all MIPs in file (1) must be sorted according to master sequence target
//coordinate, smallest to largest.
//
//This program generates genotype calls using a maximum likelihood approach with dynamic programming. For each individual, likelihoods of observing the count
//data for each MIP are calculated under 25 different possible underlying RH paralog-specific copy number states, where RHD and RHCE can have
//copy number 0-4 (25 combinations = 5*5). Then, for each individual, these likelihoods are used to construct a weighted directed acyclic graph (WDAG), with
//prior probabilities based on observed genotype data from the 1000 Genomes Project incorporated into the likelihoods for the first MIP count data. As the
//graph is constructed, highest scoring paths without any transitions between copy number states are tracked.

#include<stdio.h>
#include<string.h>
#include<gsl/gsl_randist.h>
#include<float.h>
#include<math.h>
#include<stdlib.h>
#define KRED "\x1B[31m"
#define KYEL "\x1B[33m"

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

	//get and store information about RH master sequence target coordinate for each MIP (all RH SUN-targeting MIPs have paralog-specificity)
	long start,end;
	long target_coords[num_mip_targets];
	long i=0,j;
	while(fscanf(miptargetsfile,"%s %ld %c %ld %s %c %s %s %s %s %s",dummy,&start,&miptype,&end,dummy,&miptype,dummy,dummy,dummy,dummy,dummy)==11)
  {
		if(miptype!='E')
    {
			target_coords[i]=(start+end)/2;
    	i++;
		}
  }
	fclose(miptargetsfile);

	//generate vector of possible RH copy number states and vector of ln(prior probabilities) corresponding to each copy number state	
	//assume independence of paralog-specific copy number genotypes
	int D,CE;
  int copy_states[25][2];
  double priors[25]; //vector of prior probabilities for each copy number state
		//priors_N gives the probabilities of observing a copy number state of 0,1,2,3,... for RHN
		//these probabilities are based on empirical data from 1000 Genomes Project (RH paralog-specific genotyping)
		//extremely unlikely states (never observed so far but hypothetically possible) are arbitrarily assigned a value with a resulting natural logarithm of -30
  double priors_D[5]={(double)63/(double)891,(double)274/(double)891,(double)527/(double)891,(double)27/(double)891,(exp(-30))};
  double priors_CE[5]={(exp(-30)),(double)6/(double)891,(double)866/(double)891,(double)19/(double)891,(exp(-30))};
  i=0;
	for(D=0;D<=4;D++)
  {
    for(CE=0;CE<=4;CE++)
    {
			copy_states[i][0]=D;
      copy_states[i][1]=CE;
      priors[i]=log(priors_D[D])+log(priors_CE[CE]);
      i++;
    }
  }

	struct node
	{
		long state;
		double max0;
		long index;
	};

	//allocate storage for likelihood graph
	struct node*likelihood_graph;
	likelihood_graph=(struct node*)malloc(num_mip_targets*25*sizeof(struct node));

	//initialize output files
	FILE*out;
  char output_base_name[51];
  char output_file_extension[15]=".rhsimplecalls";
	strncpy(output_base_name,*(argv+3),36);
  strcat(output_base_name,output_file_extension);
  out=fopen(output_base_name,"w");
	fprintf(out,"Individual\tRHD_CN\tRHCE_CN\tLOD_Score\n");

	//for each individual, read in counts, calculate and store individual likelihoods of data for each MIP under each copy number state
	FILE*countsfile;
	char individual[31];
	long counts[3];
	long indiv_counts[num_mip_targets][3];
	long coord,mip_coord;
	double probs[3];
	double L;
	long k;
	double mip_likelihoods[num_mip_targets][25];
	individual[30]='\0';
	countsfile=fopen(*(argv+2),"r");
	while(getc(countsfile)!='\n')
    continue;
	fgetpos(countsfile,&pos);
	while(fscanf(countsfile,"%s %s %ld %c %ld %ld %ld",individual,dummy,&coord,&miptype,&(counts[0]),&(counts[1]),&(counts[2]))==7) //while there remains data to process
	{
		fsetpos(countsfile,&pos);
		//get all relevant counts for the individual being genotyped into a single 2-dimensional array of long ([num_mip_targets][3])
		j=0;
		for(i=0;i<(num_mip_targets+num_exon_mips);i++)
		{
			fscanf(countsfile,"%s %s %ld %c %ld %ld %ld",individual,dummy,&mip_coord,&miptype,&(counts[0]),&(counts[1]),&(counts[2]));
			if(miptype!='E')
			{
				for(k=0;k<3;k++)
				{
					indiv_counts[j][k]=counts[k];
				}
				j++;
			}
		}

		//calculate likelihoods of observed counts at each MIP under each possible copy number state
		for(i=0;i<num_mip_targets;i++)
		{
			for(j=0;j<25;j++)
			{
				probs[2]=0.0;
				for(k=0;k<2;k++)
				{
					probs[k]=(double)(copy_states[j][k])/(double)(copy_states[j][0]+copy_states[j][1]);
				}
				const unsigned int countvec[3]={indiv_counts[i][0],indiv_counts[i][1],indiv_counts[i][2]};
				L=gsl_ran_multinomial_lnpdf(3,probs,countvec);
				if(L<-30)
				{
					L=-30; //minimum log likelihood value for 1 MIP probe arbitrarily assigned to -30
				}
				mip_likelihoods[i][j]=L; //[mip_target][copy_state]	
			}
		}
		
		//use dynamic programming to determine whether the data support more than 1 copy number state across SRGAP2
		double max_max0=-DBL_MAX,max2_max0=-DBL_MAX;
    long maxindex_max0;
			//setup first column and calculate maximum likelihood value and index of corresponding node
		for(j=0;j<25;j++)
    {
      likelihood_graph[j].index=j;
			likelihood_graph[j].state=j;
			likelihood_graph[j].max0=mip_likelihoods[0][j]+priors[j];
    }	
			//setup all other verticies
		for(j=25;j<(25*num_mip_targets);j++)
    {
      if((j%25)==0)
      {
        max_max0=-DBL_MAX;
        for((k=j-25);k<j;k++)
        {
          if(likelihood_graph[k].max0>max_max0)
          {
            max_max0=likelihood_graph[k].max0;
						maxindex_max0=k;
          }
        }
      }
      likelihood_graph[j].index=j;
			likelihood_graph[j].state=j%25;
			likelihood_graph[j].max0=mip_likelihoods[j/25][j%25]+likelihood_graph[j-25].max0;
    }
			//determine maximum values of final likelihoods allowing 0 transitions
		max_max0=-DBL_MAX;
		max2_max0=-DBL_MAX; //for LOD score calculation LOD=max_max0-max2_max0
    for(j=(25*num_mip_targets-25);j<(25*num_mip_targets);j++)
    {
      if(likelihood_graph[j].max0>max_max0)
      {
        if((((double)(copy_states[j%25][0])/(double)(copy_states[j%25][0]+copy_states[j%25][1]))!=((double)(copy_states[maxindex_max0][0])/(double)(copy_states[maxindex_max0][0]+copy_states[maxindex_max0][1])))||(((double)(copy_states[j%25][1])/(double)(copy_states[j%25][0]+copy_states[j%25][1]))!=((double)(copy_states[maxindex_max0][1])/(double)(copy_states[maxindex_max0][0]+copy_states[maxindex_max0][1]))))
				{
					max2_max0=max_max0; //only update value of 2nd highest likelihood if the associated copy number state has a distinct set of expected paralog-specific count frequencies from state having maximum likelihood
				}
				max_max0=likelihood_graph[j].max0;
  			maxindex_max0=j%25; //maximally likely copy number state assuming no internal events
			}
		}
		double LOD_score=max_max0-max2_max0;
		if(LOD_score>1000)
		{
			LOD_score=1000; //arbitrarily set a maximum LOD score value to 1000
		}
		fprintf(out,"%s\t%d\t%d\t%lf\n",individual,copy_states[maxindex_max0][0],copy_states[maxindex_max0][1],LOD_score);
		//prepare for next individual
		fgetpos(countsfile,&pos);
	}

	//clean up and exit
	free(likelihood_graph);
	fclose(out);
	return 0;
}